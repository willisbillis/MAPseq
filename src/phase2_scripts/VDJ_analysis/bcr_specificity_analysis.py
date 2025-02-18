import scanpy as sc
import pandas as pd
import scirpy as ir

path_data = "data"
path_bcr = f"{path_data}/BCR_01_preprocessed.h5ad"
adata_bcr = sc.read(path_bcr)

# Subset the data to include only the specified patients
patients_of_interest = ["COVID-030", "IVLPS-6", "COVID-064", "COVID-014", "COVID-027", "COVID-024"]
adata_bcr = adata_bcr[adata_bcr.obs["patient_id"].isin(patients_of_interest)].copy()

# Load the CoV-AbDab database
# The download command is commented out, assuming the file is already downloaded.
# wget -O tmp/CoV-AbDab.csv http://opig.stats.ox.ac.uk/webapps/covabdab/static/downloads/CoV-AbDab_200422.csv
cov_abdab = pd.read_csv("tmp/CoV-AbDab.csv")

# Filter for antibodies only
cov_abdab = cov_abdab[cov_abdab["Ab or Nb"] == "Ab"]

# Rename columns to match scirpy's expected format
dict_rename_cov_abdab = {
    "Heavy V Gene": "IR_VDJ_1_v_call",
    "Heavy J Gene": "IR_VDJ_1_j_call",
    "Light V Gene": "IR_VJ_1_v_call",
    "Light J Gene": "IR_VJ_1_j_call",
    "CDRH3": "IR_VDJ_1_junction_aa",
    "CDRL3": "IR_VJ_1_junction_aa",
    "Binds to": "Binding",
}
cov_abdab = cov_abdab.rename(columns=dict_rename_cov_abdab)

# Add necessary columns for compatibility with scirpy
cov_abdab["has_ir"] = "True"
cov_abdab["IR_VJ_2_junction_aa"] = None
cov_abdab["IR_VDJ_2_junction_aa"] = None

# Clean up the 'Binding' column
cov_abdab["Binding"] = cov_abdab["Binding"].apply(
    lambda x: "SARS-CoV-2" if "SARS-CoV2" in x else "None"
)

# Add cysteine and tryptophan to CDR3 sequences (important for some distance metrics)
cov_abdab["IR_VJ_1_junction_aa"] = "C" + cov_abdab["IR_VJ_1_junction_aa"] + "F"
cov_abdab["IR_VDJ_1_junction_aa"] = "C" + cov_abdab["IR_VDJ_1_junction_aa"] + "W"

# Convert the database to an AnnData object
cov_abdab = sc.AnnData(obs=cov_abdab)
cov_abdab.uns["DB"] = {}
cov_abdab.uns["DB"]["name"] = "CoV-AbDab"

# Perform matching using different metrics

# Identity matching
metric = "identity"
sequence = "aa"

ir.pp.ir_dist(adata_bcr, cov_abdab, metric=metric, sequence=sequence)
ir.tl.ir_query(
    adata_bcr,
    cov_abdab,
    metric=metric,
    sequence=sequence,
    receptor_arms="VDJ",
    dual_ir="primary_only",
)
ir.tl.ir_query_annotate(
    adata_bcr, cov_abdab, metric=metric, sequence=sequence, include_ref_cols=["Binding"]
)
print(adata_bcr.obs["Binding"].value_counts())


# Identity matching with all receptor arms
ir.tl.ir_query(
    adata_bcr,
    cov_abdab,
    metric=metric,
    sequence=sequence,
    receptor_arms="all",
    dual_ir="primary_only",
)
ir.tl.ir_query_annotate(
    adata_bcr,
    cov_abdab,
    metric=metric,
    sequence=sequence,
    include_ref_cols=["Binding"],
    suffix="_fullIR",
)
print(adata_bcr.obs["Binding_fullIR"].value_counts())

# Hamming distance matching
metric = "hamming"
sequence = "aa"

ir.pp.ir_dist(adata_bcr, cov_abdab, metric=metric, sequence=sequence)
ir.tl.ir_query(
    adata_bcr,
    cov_abdab,
    metric=metric,
    sequence=sequence,
    receptor_arms="all",
    dual_ir="primary_only",
)
ir.tl.ir_query_annotate(
    adata_bcr,
    cov_abdab,
    metric=metric,
    sequence=sequence,
    include_ref_cols=["Binding"],
    suffix="_Hamming",
)
print(adata_bcr.obs["Binding_Hamming"].value_counts())

