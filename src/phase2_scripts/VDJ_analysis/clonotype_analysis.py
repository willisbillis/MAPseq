import warnings

warnings.filterwarnings(
    "ignore",
    ".*IProgress not found*",
)
warnings.simplefilter(action="ignore", category=FutureWarning)

import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir
from palmotif import compute_motif, svg_logo

warnings.simplefilter(action="ignore", category=pd.errors.DtypeWarning)

path_data = "/home/icb/juan.henao/BestPracticeStart/data"

path_gex = f"{path_data}/TCR_filtered.h5ad"
adata = sc.read(path_gex)

adata = adata[
    adata.obs["patient_id"].isin(
        ["COVID-014", "CV0902", "AP6", "COVID-045", "COVID-066", "COVID-067"]
    )
]
# adata[adata.obs['Status'] == 'Healthy'].obs['patient_id']
_ = ir.pl.group_abundance(
    adata,
    groupby="patient_id",
    target_col="chain_pairing",
    normalize=True,
    figsize=[10, 10],
)

ir.pp.ir_dist(adata, sequence="aa")
ir.tl.define_clonotype_clusters(
    adata, sequence="aa", receptor_arms="all", dual_ir="primary_only"
)
ir.tl.clonotype_network(adata, min_cells=50, sequence="aa")
_ = ir.pl.clonotype_network(
    adata,
    color="patient_id",
    base_size=10,
    label_fontsize=9,
    panel_size=(10, 10),
    legend_fontsize=15,
)

adata.obs["cc_aa_identity"] = adata.obs["cc_aa_identity"].astype("str")
adata.obs.loc[adata.obs["cc_aa_identity"] == "0", :].groupby(
    [
        "IR_VJ_1_junction_aa",
        "IR_VDJ_1_junction_aa",
        "receptor_subtype",
    ],
    observed=True,
).size().reset_index(name="n_cells")

ir.tl.clonal_expansion(adata, target_col="cc_aa_identity")

_ = ir.pl.clonal_expansion(
    adata,
    groupby="initial_clustering",
    target_col="cc_aa_identity",
    clip_at=4,
    normalize=False,
    figsize=[10, 10],
)
_ = ir.pl.clonal_expansion(
    adata, "initial_clustering", target_col="cc_aa_identity", figsize=[10, 10]
)
_ = ir.pl.alpha_diversity(
    adata, groupby="initial_clustering", target_col="cc_aa_identity", figsize=[10, 10]
)
_ = ir.pl.group_abundance(
    adata,
    groupby="cc_aa_identity",
    target_col="initial_clustering",
    max_cols=10,
    figsize=[10, 10],
)
# By condition

_ = ir.pl.group_abundance(
    adata,
    groupby="cc_aa_identity",
    target_col="Status",
    max_cols=15,
    fig_kws={"dpi": 100},
    figsize=[10, 10],
)

# By sample

_ = ir.pl.group_abundance(
    adata,
    groupby="cc_aa_identity",
    target_col="patient_id",
    max_cols=15,
    fig_kws={"dpi": 100},
    figsize=[10, 10],
)

_ = ir.pl.group_abundance(
    adata,
    groupby="IR_VJ_1_v_call",
    target_col="initial_clustering",
    normalize=False,
    max_cols=10,
    figsize=[10, 10],
)

# Normalized abundances

_ = ir.pl.group_abundance(
    adata,
    groupby="IR_VJ_1_v_call",
    target_col="initial_clustering",
    normalize=True,
    max_cols=10,
    figsize=[10, 10],
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
)

_ = ir.pl.vdj_usage(
    adata,
    full_combination=False,
    max_segments=None,
    max_ribbons=30,
    fig_kws={"figsize": [10, 10]},
)

adata.obs[adata.obs["IR_VDJ_1_d_call"] == "TRBD2"].cc_aa_identity.value_counts()

_ = ir.pl.vdj_usage(
    adata[
        adata.obs["cc_aa_identity"].isin(["5150", "1815", "4427", "3500", "4078"]), :
    ],
    max_ribbons=None,
    max_segments=100,
    fig_kws={"figsize": [10, 10]},
)

_ = ir.pl.spectratype(
    adata,
    cdr3_col="IR_VDJ_1_junction_aa",
    color="initial_clustering",
    viztype="bar",
    fig_kws={"dpi": 120},
    figsize=[10, 10],
)

_ = ir.pl.spectratype(
    adata,
    cdr3_col="IR_VDJ_1_junction_aa",
    color="initial_clustering",
    viztype="curve",
    curve_layout="shifted",
    fig_kws={"figsize": [10, 10]},
    kde_kws={"kde_norm": False},
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
)

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

df, dst, lk = ir.tl.repertoire_overlap(
    adata, "patient_id", target_col="cc_aa_identity", inplace=False
)

ir.pl.repertoire_overlap(
    adata, "patient_id", target_col="cc_aa_identity", heatmap_cats=["Centre"]
)

_ = ir.pl.repertoire_overlap(
    adata,
    "patient_id",
    pair_to_plot=["COVID-067", "CV0902"],
    fig_kws={"figsize": [10, 10]},
)