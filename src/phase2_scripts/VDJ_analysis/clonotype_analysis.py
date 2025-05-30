import scanpy as sc
import scirpy as ir
from palmotif import compute_motif, svg_logo

# Define data paths
path_data = "/home/icb/juan.henao/BestPracticeStart/data"
path_gex = f"{path_data}/TCR_filtered.h5ad"

# Load data. Handle potential exceptions (e.g., FileNotFoundError).
try:
    adata = sc.read(path_gex)
except FileNotFoundError:
    print(f"Error: File not found at {path_gex}")
    exit(1)


# Subset data to specific patients.  Consider making this configurable.
adata = adata[
    adata.obs["patient_id"].isin(
        ["COVID-014", "CV0902", "AP6", "COVID-045", "COVID-066", "COVID-067"]
    )
]

# Analyze and visualize clonotypes

#Visualize the abundance of different chain pairings per patient
_ = ir.pl.group_abundance(
    adata,
    groupby="patient_id",
    target_col="chain_pairing",
    normalize=True,
    figsize=[10, 10],
    title="Chain Pairing Abundance per Patient"
)

#Compute the distance between Immunoglobulin/T-cell receptor sequences
ir.pp.ir_dist(adata, sequence="aa")

#Define clonotype clusters based on amino acid sequence
ir.tl.define_clonotype_clusters(
    adata, sequence="aa", receptor_arms="all", dual_ir="primary_only"
)

#Build a clonotype network
ir.tl.clonotype_network(adata, min_cells=50, sequence="aa")

#Visualize the clonotype network colored by patient ID.
_ = ir.pl.clonotype_network(
    adata,
    color="patient_id",
    base_size=10,
    label_fontsize=9,
    panel_size=(10, 10),
    legend_fontsize=15,
    title="Clonotype Network by Patient ID"
)

#Data type conversion and analysis of clonotypes with no chains
adata.obs["cc_aa_identity"] = adata.obs["cc_aa_identity"].astype("str")
clonotype_no_chains = adata.obs.loc[adata.obs["cc_aa_identity"] == "0", :].groupby(
    [
        "IR_VJ_1_junction_aa",
        "IR_VDJ_1_junction_aa",
        "receptor_subtype",
    ],
    observed=True,
).size().reset_index(name="n_cells")
print("Clonotypes with no chains:\n", clonotype_no_chains) #Print the result for better understanding


# Analyze clonal expansion
ir.tl.clonal_expansion(adata, target_col="cc_aa_identity")

#Visualize clonal expansion by initial clustering
_ = ir.pl.clonal_expansion(
    adata,
    groupby="initial_clustering",
    target_col="cc_aa_identity",
    clip_at=4,
    normalize=False,
    figsize=[10, 10],
    title="Clonal Expansion by Initial Clustering"
)

_ = ir.pl.clonal_expansion(
    adata, "initial_clustering", target_col="cc_aa_identity", figsize=[10, 10], title="Clonal Expansion"
)

_ = ir.pl.alpha_diversity(
    adata, groupby="initial_clustering", target_col="cc_aa_identity", figsize=[10, 10], title="Alpha Diversity"
)

_ = ir.pl.group_abundance(
    adata,
    groupby="cc_aa_identity",
    target_col="initial_clustering",
    max_cols=10,
    figsize=[10, 10],
    title="Clonotype Abundance by Initial Clustering"
)


# Group abundance plots by condition and sample.
_ = ir.pl.group_abundance(
    adata,
    groupby="cc_aa_identity",
    target_col="Status",
    max_cols=15,
    fig_kws={"dpi": 100},
    figsize=[10, 10],
    title="Clonotype Abundance by Status"
)

_ = ir.pl.group_abundance(
    adata,
    groupby="cc_aa_identity",
    target_col="patient_id",
    max_cols=15,
    fig_kws={"dpi": 100},
    figsize=[10, 10],
    title="Clonotype Abundance by Patient ID"
)

_ = ir.pl.group_abundance(
    adata,
    groupby="IR_VJ_1_v_call",
    target_col="initial_clustering",
    normalize=False,
    max_cols=10,
    figsize=[10, 10],
    title="V-gene Usage by Initial Clustering (Unnormalized)"
)

_ = ir.pl.group_abundance(
    adata,
    groupby="IR_VJ_1_v_call",
    target_col="initial_clustering",
    normalize=True,
    max_cols=10,
    figsize=[10, 10],
    title="V-gene Usage by Initial Clustering (Normalized)"
)

_ = ir.pl.group_abundance(
    adata[
        adata.obs["IR_VDJ_1_v_call"].isin(
            ["TRBV19", "TRBV10-1", "TRBV11-1", "TRBV7-9"]
        ),
        :,
    ],
    groupby="initial_clustering",
    target_col="IR_VDJ_1_v_call",
    normalize=True,
    figsize=[10, 10],
    title="Specific V-gene Usage by Initial Clustering"
)

_ = ir.pl.vdj_usage(
    adata,
    full_combination=False,
    max_segments=None,
    max_ribbons=30,
    fig_kws={"figsize": [10, 10]},
    title="VDJ Usage"
)

print("Count of cc_aa_identity where IR_VDJ_1_d_call == 'TRBD2':\n", adata.obs[adata.obs["IR_VDJ_1_d_call"] == "TRBD2"].cc_aa_identity.value_counts())

_ = ir.pl.vdj_usage(
    adata[
        adata.obs["cc_aa_identity"].isin(["5150", "1815", "4427", "3500", "4078"]), :
    ],
    max_ribbons=None,
    max_segments=100,
    fig_kws={"figsize": [10, 10]},
    title="VDJ Usage for Specific Clonotypes"
)


_ = ir.pl.spectratype(
    adata,
    cdr3_col="IR_VDJ_1_junction_aa",
    color="initial_clustering",
    viztype="bar",
    fig_kws={"dpi": 120},
    figsize=[10, 10],
    title="Spectratype (Bar)"
)

_ = ir.pl.spectratype(
    adata,
    cdr3_col="IR_VDJ_1_junction_aa",
    color="initial_clustering",
    viztype="curve",
    curve_layout="shifted",
    fig_kws={"figsize": [10, 10]},
    kde_kws={"kde_norm": False},
    title="Spectratype (Curve)"
)

_ = ir.pl.spectratype(
    adata[
        adata.obs["IR_VDJ_1_v_call"].isin(
            ["TRBV5-1", "TRBV11-2", "TRBV7-2", "TRBV11-3"]
        ),
        :,
    ],
    cdr3_col="IR_VDJ_1_junction_aa",
    color="IR_VDJ_1_v_call",
    normalize="initial_clustering",
    fig_kws={"dpi": 120},
    figsize=[10, 10],
    title="Spectratype for Specific V-genes"
)

#Compute the motif. Handle potential errors (e.g., empty input)
try:
    motif = compute_motif(
        adata[
            (
                adata.obs["IR_VDJ_1_v_call"].isin(
                    ["TRBV5-1", "TRBV11-2", "TRBV7-2", "TRBV11-3"]
                )
            )
            & (adata.obs["IR_VDJ_1_junction_aa"].str.len() == 15),
            :,
        ]
        .obs["IR_VDJ_1_junction_aa"]
        .to_list()
    )
    _ = svg_logo(
        motif, "../_static/images/air_repertoire/logo_motif.svg", color_scheme="taylor"
    )
except ValueError as e:
    print(f"Error computing motif: {e}")


# Analyze repertoire overlap
df, dst, lk = ir.tl.repertoire_overlap(
    adata, "patient_id", target_col="cc_aa_identity", inplace=False
)

ir.pl.repertoire_overlap(
    adata, "patient_id", target_col="cc_aa_identity", heatmap_cats=["Centre"], title="Repertoire Overlap"
)

_ = ir.pl.repertoire_overlap(
    adata,
    "patient_id",
    pair_to_plot=["COVID-067", "CV0902"],
    fig_kws={"figsize": [10, 10]},
    title="Repertoire Overlap (COVID-067 vs CV0902)"
)

