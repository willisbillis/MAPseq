import warnings

warnings.filterwarnings(
    "ignore",
    ".*IProgress not found*",
)

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scirpy as ir
import seaborn as sb

path_data = "data"
path_tcr = f"{path_data}/TCR_01_preprocessed.h5ad"
adata_tcr = sc.read(path_tcr)

# TODO final decision sampling: prosponed until remaining chapters are fixed
adata_tcr = adata_tcr[adata_tcr.obs["patient_id"].isin(["CV0902", "AP6"])].copy()

vdjdb = ir.datasets.vdjdb()
vdjdb.obs.head(5)

vdjdb[vdjdb.obs["antigen.species"] == "SARS-CoV-2"].obs[
    "antigen.epitope"
].value_counts()

vdjdb[vdjdb.obs["antigen.species"] == "SARS-CoV-2"].obs[
    "antigen.species"
].value_counts()

# Manual Query
print(f"Amount of samples in VDJDB: {len(vdjdb)}")
vdjdb_overlap = vdjdb[
    vdjdb.obs["IR_VDJ_1_junction_aa"].isin(adata_tcr.obs["IR_VDJ_1_junction_aa"])
].obs
vdjdb_overlap = vdjdb_overlap[["IR_VDJ_1_junction_aa", "antigen.species"]]
print(f"Amount of overlapping samples in VDJDB: {len(vdjdb_overlap)}")

adata_tcr.obs["has_vdjdb_overlap"] = adata_tcr.obs["IR_VDJ_1_junction_aa"].isin(
    vdjdb_overlap["IR_VDJ_1_junction_aa"]
)
adata_tcr.obs["has_vdjdb_overlap"].value_counts()

def assign_disease(cdr3beta):
    matching_rows = vdjdb_overlap[vdjdb_overlap["IR_VDJ_1_junction_aa"] == cdr3beta]
    diseases = matching_rows["antigen.species"].values.tolist()
    diseases = list(set(diseases))
    if len(diseases) > 1:
        return "ambiguous"
    if len(diseases) == 0:
        return "no entry"
    return diseases[0]

adata_tcr.obs["antigen.species_manual"] = (
    adata_tcr[adata_tcr.obs["has_vdjdb_overlap"]]
    .obs["IR_VDJ_1_junction_aa"]
    .apply(assign_disease)
)
adata_tcr.obs["antigen.species_manual"].value_counts()

# Query with various strictness
metric = "identity"
sequence = "aa"
ir.pp.ir_dist(adata_tcr, vdjdb, metric=metric, sequence=sequence)

ir.tl.ir_query(
    adata_tcr,
    vdjdb,
    metric=metric,
    sequence=sequence,
    receptor_arms="VDJ",
    dual_ir="primary_only",
)

ir.tl.ir_query_annotate(
    adata_tcr,
    vdjdb,
    metric=metric,
    sequence=sequence,
    include_ref_cols=["antigen.species", "antigen.epitope"],
)
adata_tcr.obs["antigen.species"].value_counts()

ir.tl.ir_query(
    adata_tcr,
    vdjdb,
    metric=metric,
    sequence=sequence,
    receptor_arms="VJ",
    dual_ir="primary_only",
)
ir.tl.ir_query_annotate(
    adata_tcr,
    vdjdb,
    metric=metric,
    sequence=sequence,
    include_ref_cols=["antigen.species", "antigen.epitope"],
    suffix="_VJ",
)
adata_tcr.obs["antigen.species_VJ"].value_counts()

ir.tl.ir_query(
    adata_tcr,
    vdjdb,
    metric=metric,
    sequence=sequence,
    receptor_arms="all",
    dual_ir="primary_only",
)
ir.tl.ir_query_annotate(
    adata_tcr,
    vdjdb,
    metric=metric,
    sequence=sequence,
    include_ref_cols=["antigen.species", "antigen.epitope"],
    suffix="_fullIR",
)
adata_tcr.obs["antigen.species_fullIR"].value_counts()

## TCRdist

from tcrdist.repertoire import TCRrep

adata_tcrdist = adata_tcr[
    ~adata_tcr.obs["chain_pairing"].isin(["No IR", "orphan VDJ", "orphan VJ"])
]

for col in ["IR_VJ_1_v_call", "IR_VDJ_1_v_call"]:
    adata_tcrdist = adata_tcrdist[~adata_tcrdist.obs[col].isna()]
adata_tcrdist = adata_tcrdist.copy()


df_tcrdist = adata_tcrdist.obs[
    [
        "IR_VJ_1_junction_aa",
        "IR_VJ_1_v_call",
        "IR_VJ_1_j_call",
        "IR_VDJ_1_junction_aa",
        "IR_VDJ_1_v_call",
        "IR_VDJ_1_j_call",
        "chain_pairing",
        "patient_id",
        "antigen.species",
        "initial_clustering",
    ]
].copy()

for col in df_tcrdist.columns:
    df_tcrdist[col] = df_tcrdist[col].astype(str)
df_tcrdist = df_tcrdist.drop_duplicates()
df_tcrdist = df_tcrdist.reset_index(drop=True)

dict_rename_tcrdist = {
    "IR_VJ_1_junction_aa": "cdr3_a_aa",
    "IR_VJ_1_v_call": "v_a_gene",
    "IR_VDJ_1_junction_aa": "cdr3_b_aa",
    "IR_VDJ_1_v_call": "v_b_gene",
}

df_tcrdist = df_tcrdist.rename(columns=dict_rename_tcrdist)

df_tcrdist["v_a_gene"] = df_tcrdist["v_a_gene"].astype(str) + "*01"
df_tcrdist["v_b_gene"] = df_tcrdist["v_b_gene"].astype(str) + "*01"

df_tcrdist["count"] = 1  # todo

df_tcrdist["index"] = df_tcrdist.index

df_tcrdist.head(5)

tr = TCRrep(cell_df=df_tcrdist, organism="human", chains=["alpha", "beta"])

dist_total = tr.pw_alpha + tr.pw_beta
columns = tr.clone_df["index"].astype(float).astype(int)
df_dist = pd.DataFrame(dist_total, columns=columns, index=columns)

import scipy.cluster.hierarchy as hc
import scipy.spatial as sp
from matplotlib import rcParams

rcParams["figure.figsize"] = (20, 20)

linkage = hc.linkage(sp.distance.squareform(df_dist), method="average")
plot = sb.clustermap(
    df_dist,
    row_linkage=linkage,
    col_linkage=linkage,
    figsize=(40, 40),
    yticklabels=False,
    xticklabels=False,
)
plot.ax_row_dendrogram.set_visible(False)
plot.ax_col_dendrogram.set_visible(False)

plot.figure.suptitle("TCRdist")
plt.tight_layout()

df_tcrdist_alpha = pd.DataFrame(
    tr.pw_alpha, columns=tr.clone_df["cdr3_a_aa"], index=tr.clone_df["cdr3_a_aa"]
)
df_tcrdist_beta = pd.DataFrame(
    tr.pw_beta, columns=tr.clone_df["cdr3_b_aa"], index=tr.clone_df["cdr3_b_aa"]
)

from scipy.sparse import csr_matrix


def add_dists(adata, df_dist_alpha, df_dist_beta, name, cutoff):
    adata.uns[f"ir_dist_aa_{name}"] = {
        "params": {"metric": f"{name}", "sequence": "aa", "cutoff": cutoff}
    }

    for chain, dists in [("VJ", df_dist_alpha), ("VDJ", df_dist_beta)]:
        if dists is None:
            continue
        dist_values = dists.values + 1
        dist_values[dist_values > (cutoff + 1)] = 0
        dist_values = csr_matrix(dist_values)
        adata.uns[f"ir_dist_aa_{name}"][chain] = {
            "seqs": dists.index.tolist(),
            "distances": dist_values,
        }


add_dists(adata_tcrdist, df_tcrdist_alpha, df_tcrdist_beta, "tcrdist", 60)

ir.tl.define_clonotype_clusters(
    adata_tcrdist,
    sequence="aa",
    metric="tcrdist",
    receptor_arms="all",
    dual_ir="primary_only",
)

adata_tcrdist.obs["antigen.species"] = adata_tcrdist.obs["antigen.species"].astype(str)

ir.tl.clonotype_network(
    adata_tcrdist, min_cells=3, min_nodes=3, sequence="aa", metric="tcrdist"
)
ir.pl.clonotype_network(
    adata_tcrdist,
    color="antigen.species",
    label_fontsize=9,
    panel_size=(14, 7),
    base_size=10,
    size_power=0.75,
)