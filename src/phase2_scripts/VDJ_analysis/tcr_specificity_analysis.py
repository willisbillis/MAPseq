import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scirpy as ir
import seaborn as sb
import scipy.cluster.hierarchy as hc
import scipy.spatial as sp
from matplotlib import rcParams
from scipy.sparse import csr_matrix
from tcrdist.repertoire import TCRrep


# Define paths
path_data = "data"
path_tcr = f"{path_data}/TCR_01_preprocessed.h5ad"

# Load TCR data
adata_tcr = sc.read(path_tcr)

# Subset data for specific patients (optional, remove if not needed)
adata_tcr = adata_tcr[adata_tcr.obs["patient_id"].isin(["CV0902", "AP6"])].copy()

# Load VDJdb data
vdjdb = ir.datasets.vdjdb()

# Explore VDJdb data (optional)
# print(f"Amount of samples in VDJDB: {len(vdjdb)}")
# print(vdjdb[vdjdb.obs['antigen.species'] == 'SARS-CoV-2'].obs['antigen.epitope'].value_counts())
# print(vdjdb[vdjdb.obs['antigen.species'] == 'SARS-CoV-2'].obs['antigen.species'].value_counts())


# Manual VDJdb overlap analysis (less robust than scirpy's methods, consider removing)
vdjdb_overlap = vdjdb[vdjdb.obs["IR_VDJ_1_junction_aa"].isin(adata_tcr.obs["IR_VDJ_1_junction_aa"])].obs
vdjdb_overlap = vdjdb_overlap[["IR_VDJ_1_junction_aa", "antigen.species"]]
# print(f"Amount of overlapping samples in VDJDB: {len(vdjdb_overlap)}")

adata_tcr.obs["has_vdjdb_overlap"] = adata_tcr.obs["IR_VDJ_1_junction_aa"].isin(vdjdb_overlap["IR_VDJ_1_junction_aa"])


def assign_disease(cdr3beta):
    matching_rows = vdjdb_overlap[vdjdb_overlap["IR_VDJ_1_junction_aa"] == cdr3beta]
    diseases = matching_rows["antigen.species"].unique().tolist() # Use unique to avoid duplicates
    if len(diseases) > 1:
        return "ambiguous"
    if len(diseases) == 0:
        return "no entry"
    return diseases[0]


adata_tcr.obs["antigen.species_manual"] = adata_tcr[adata_tcr.obs["has_vdjdb_overlap"]].obs["IR_VDJ_1_junction_aa"].apply(assign_disease)


# Scirpy-based VDJdb queries (recommended approach)
metric = "identity"
sequence = "aa"

# VDJ query
ir.pp.ir_dist(adata_tcr, vdjdb, metric=metric, sequence=sequence)
ir.tl.ir_query(adata_tcr, vdjdb, metric=metric, sequence=sequence, receptor_arms="VDJ", dual_ir="primary_only")
ir.tl.ir_query_annotate(adata_tcr, vdjdb, metric=metric, sequence=sequence, include_ref_cols=["antigen.species", "antigen.epitope"])

# VJ query (optional)
ir.tl.ir_query(adata_tcr, vdjdb, metric=metric, sequence=sequence, receptor_arms="VJ", dual_ir="primary_only")
ir.tl.ir_query_annotate(adata_tcr, vdjdb, metric=metric, sequence=sequence, include_ref_cols=["antigen.species", "antigen.epitope"], suffix="_VJ")

# All chains query (optional)
ir.tl.ir_query(adata_tcr, vdjdb, metric=metric, sequence=sequence, receptor_arms="all", dual_ir="primary_only")
ir.tl.ir_query_annotate(adata_tcr, vdjdb, metric=metric, sequence=sequence, include_ref_cols=["antigen.species", "antigen.epitope"], suffix="_fullIR")


# TCRdist analysis
adata_tcrdist = adata_tcr[~adata_tcr.obs["chain_pairing"].isin(["No IR", "orphan VDJ", "orphan VJ"])].copy()  # Create a copy to avoid modifying the original adata object

# Filter out cells with missing V/J gene information
for col in ["IR_VJ_1_v_call", "IR_VDJ_1_v_call"]:
    adata_tcrdist = adata_tcrdist[~adata_tcrdist.obs[col].isna()]


# Prepare data for TCRdist
df_tcrdist = adata_tcrdist.obs[
    [
        "IR_VJ_1_junction_aa", "IR_VJ_1_v_call", "IR_VDJ_1_junction_aa",
        "IR_VDJ_1_v_call", "chain_pairing", "patient_id", "antigen.species", "initial_clustering"
    ]
].copy()

# Convert columns to string type
for col in df_tcrdist.columns:
    df_tcrdist[col] = df_tcrdist[col].astype(str)

# Drop duplicate TCR sequences
df_tcrdist = df_tcrdist.drop_duplicates().reset_index(drop=True)


# Rename columns to match TCRrep format
dict_rename_tcrdist = {
    "IR_VJ_1_junction_aa": "cdr3_a_aa", "IR_VJ_1_v_call": "v_a_gene",
    "IR_VDJ_1_junction_aa": "cdr3_b_aa", "IR_VDJ_1_v_call": "v_b_gene",
}
df_tcrdist = df_tcrdist.rename(columns=dict_rename_tcrdist)


# Add allele information (required by TCRrep, using *01 as placeholder)
df_tcrdist["v_a_gene"] = df_tcrdist["v_a_gene"].astype(str) + "*01"
df_tcrdist["v_b_gene"] = df_tcrdist["v_b_gene"].astype(str) + "*01"


# Add count column (required by TCRrep, assuming 1 for each unique sequence)
df_tcrdist["count"] = 1

# Add index column
df_tcrdist["index"] = df_tcrdist.index


# Create TCRrep object
tr = TCRrep(cell_df=df_tcrdist, organism="human", chains=["alpha", "beta"])


# Calculate TCR distances
dist_total = tr.pw_alpha + tr.pw_beta
columns = tr.clone_df["index"].astype(float).astype(int)
df_dist = pd.DataFrame(dist_total, columns=columns, index=columns)


# Clustermap visualization (optional)
rcParams["figure.figsize"] = (20, 20)
linkage = hc.linkage(sp.distance.squareform(df_dist), method="average")

plot = sb.clustermap(df_dist, row_linkage=linkage, col_linkage=linkage, figsize=(40, 40), yticklabels=False, xticklabels=False)
plot.ax_row_dendrogram.set_visible(False)
plot.ax_col_dendrogram.set_visible(False)
plot.figure.suptitle("TCRdist")
plt.tight_layout()


# Calculate distance matrices for alpha and beta chains
df_tcrdist_alpha = pd.DataFrame(tr.pw_alpha, columns=tr.clone_df["cdr3_a_aa"], index=tr.clone_df["cdr3_a_aa"])
df_tcrdist_beta = pd.DataFrame(tr.pw_beta, columns=tr.clone_df["cdr3_b_aa"], index=tr.clone_df["cdr3_b_aa"])


# Function to add TCRdist distances to AnnData object
def add_dists(adata, df_dist_alpha, df_dist_beta, name, cutoff):
    adata.uns[f"ir_dist_aa_{name}"] = {"params": {"metric": f"{name}", "sequence": "aa", "cutoff": cutoff}}
    for chain, dists in [("VJ", df_dist_alpha), ("VDJ", df_dist_beta)]:
        if dists is None:
            continue
        dist_values = dists.values + 1
        dist_values[dist_values > (cutoff + 1)] = 0  # Apply cutoff
        dist_values = csr_matrix(dist_values)  # Convert to sparse matrix
        adata.uns[f"ir_dist_aa_{name}"][chain] = {"seqs": dists.index.tolist(), "distances": dist_values}


# Add TCRdist distances to AnnData object
add_dists(adata_tcrdist, df_tcrdist_alpha, df_tcrdist_beta, "tcrdist", 60)


# Define clonotype clusters based on TCRdist
ir.tl.define_clonotype_clusters(adata_tcrdist, sequence="aa", metric="tcrdist", receptor_arms="all", dual_ir="primary_only")

adata_tcrdist.obs["antigen.species"] = adata_tcrdist.obs["antigen.species"].astype(str)


# Create and plot clonotype network
ir.tl.clonotype_network(adata_tcrdist, min_cells=3, min_nodes=3, sequence="aa", metric="tcrdist")
ir.pl.clonotype_network(adata_tcrdist, color="antigen.species", label_fontsize=9, panel_size=(14, 7), base_size=10, size_power=0.75)

# Save the updated anndata object (optional)
# adata_tcrdist.write_h5ad("updated_tcr_data.h5ad")

plt.show() # Necessary to display plots when running in a script
