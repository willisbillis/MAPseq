import scanpy as sc
import pandas as pd
import scirpy as ir

path_data = "data"
path_bcr = f"{path_data}/BCR_01_preprocessed.h5ad"
adata_bcr = sc.read(path_bcr)

# Subset the data to include only the specified patients
patients_of_interest = ["COVID-030", "IVLPS-6", "COVID-064", "COVID-014", "COVID-027", "COVID-024"]
adata_bcr = adata_bcr[adata_bcr.obs["patient_id"].isin(patients_of_interest)].copy()

# Load the influenza antibody database
# The download command is commented out, assuming the file is already downloaded.
# wget -O tmp/Influenza-AbDab.csv http://example.com/influenza_abdab.csv
influenza_abdab = pd.read_csv("tmp/Influenza-AbDab.csv")

# Filter for antibodies only
influenza_abdab = influenza_abdab[influenza_abdab["Ab or Nb"] == "Ab"]

# Rename columns to match scirpy's expected format
dict_rename_influenza_abdab = {
    "Heavy V Gene": "IR_VDJ_1_v_call",
    "Heavy J Gene": "IR_VDJ_1_j_call",
    "Light V Gene": "IR_VJ_1_v_call",
    "Light J Gene": "IR_VJ_1_j_call",
    "CDRH3": "IR_VDJ_1_junction_aa",
    "CDRL3": "IR_VJ_1_junction_aa",
    "Binds to": "Binding",
}
influenza_abdab = influenza_abdab.rename(columns=dict_rename_influenza_abdab)

# Add necessary columns for compatibility with scirpy
influenza_abdab["has_ir"] = "True"
influenza_abdab["IR_VJ_2_junction_aa"] = None
influenza_abdab["IR_VDJ_2_junction_aa"] = None

# Clean up the 'Binding' column
influenza_abdab["Binding"] = influenza_abdab["Binding"].apply(
    lambda x: "Influenza" if isinstance(x, str) and "Influenza" in x else ("None" if pd.isna(x) or not isinstance(x, str) else x)
)

# Add cysteine and tryptophan to CDR3 sequences (important for some distance metrics)
# Fill missing CDR3 sequences with empty strings before concatenation
influenza_abdab["IR_VJ_1_junction_aa"] = influenza_abdab["IR_VJ_1_junction_aa"].fillna("")
influenza_abdab["IR_VDJ_1_junction_aa"] = influenza_abdab["IR_VDJ_1_junction_aa"].fillna("")
influenza_abdab["IR_VJ_1_junction_aa"] = "C" + influenza_abdab["IR_VJ_1_junction_aa"].astype(str) + "F"
influenza_abdab["IR_VDJ_1_junction_aa"] = "C" + influenza_abdab["IR_VDJ_1_junction_aa"].astype(str) + "W"
influenza_abdab.uns["DB"] = {"name": "Influenza-AbDab"}

# Perform matching using different metrics

# Identity matching
metric = "identity"
sequence = "aa"

ir.pp.ir_dist(adata_bcr, influenza_abdab, metric=metric, sequence=sequence)
ir.tl.ir_query(
    adata_bcr,
    influenza_abdab,
    metric=metric,
    sequence=sequence,
    receptor_arms="VDJ",
    dual_ir="primary_only",
)
ir.tl.ir_query_annotate(
    adata_bcr, influenza_abdab, metric=metric, sequence=sequence, include_ref_cols=["Binding"]
)
print(adata_bcr.obs["Binding"].value_counts())


# Identity matching with all receptor arms
ir.tl.ir_query(
    adata_bcr,
    influenza_abdab,
    metric=metric,
    sequence=sequence,
    receptor_arms="all",
    dual_ir="primary_only",
)
ir.tl.ir_query_annotate(
    adata_bcr,
    influenza_abdab,
    metric=metric,
    sequence=sequence,
    include_ref_cols=["Binding"],
    suffix="_fullIR",
)
print(adata_bcr.obs["Binding_fullIR"].value_counts())

# Hamming distance matching
metric = "hamming"
sequence = "aa"

ir.pp.ir_dist(adata_bcr, influenza_abdab, metric=metric, sequence=sequence)
ir.tl.ir_query(
    adata_bcr,
    influenza_abdab,
    metric=metric,
    sequence=sequence,
    receptor_arms="all",
    dual_ir="primary_only",
)
ir.tl.ir_query_annotate(
    adata_bcr,
    influenza_abdab,
    metric=metric,
    sequence=sequence,
    include_ref_cols=["Binding"],
    suffix="_Hamming",
)
print(adata_bcr.obs["Binding_Hamming"].value_counts())
