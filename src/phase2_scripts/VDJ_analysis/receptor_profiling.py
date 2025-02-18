import warnings

warnings.filterwarnings(
    "ignore",
    ".*IProgress not found*",
)
warnings.simplefilter(action="ignore", category=FutureWarning)

import pandas as pd
import scanpy as sc
import scirpy as ir

warnings.simplefilter(action="ignore", category=pd.errors.DtypeWarning)

path_data = "data/"
path_bcr_input = f"{path_data}/BCR_00_read_aligned.csv"
path_tcr_input = f"{path_data}/TCR_00_read_aligned.tsv"

path_bcr_out = f"{path_data}/BCR_01_preprocessed.h5ad"

path_tcr_csv = f"{path_data}/TCR_00_read_aligned.csv"
path_tcr_out = f"{path_data}/TCR_01_preprocessed.h5ad"

df_bcr_raw = pd.read_csv(path_bcr_input, index_col=0)

# The column 'productive' contains mixed data types which are not compatible with downstream tools.
# We correct them by casting them to strings.
df_bcr_raw["productive"] = df_bcr_raw["productive"].astype(str)
print(f"Total measurements: {len(df_bcr_raw)}")
df_bcr_raw.head(5)

columns = [
    "barcode",
    "v_gene",
    "d_gene",
    "c_gene",
    "j_gene",
    "productive",
    "full_length",
]

columns += ["cdr3", "cdr3_nt"]

df_tcr_raw = pd.read_csv(path_tcr_input, sep="\t")
df_tcr_raw["barcode"] = df_tcr_raw.pop("CellID")
df_tcr_raw.to_csv(path_tcr_csv)
print(f"Total measurements: {len(df_tcr_raw)}")

adata_tcr = ir.io.read_10x_vdj(path_tcr_csv)
print(f"Amount cells: {len(adata_tcr)}")

patient_information = [
    "barcode",
    "Centre",
    "Sample",
    "patient_id",
    "Collection_Day",
    "Sex",
    "Swab_result",
    "Status",
    "Smoker",
    "Status_on_day_collection",
    "Status_on_day_collection_summary",
    "Days_from_onset",
    "time_after_LPS",
    "Worst_Clinical_Status",
    "Outcome",
    "initial_clustering",
    "study_id",
    "AgeRange",
    "Age",
]
df_patient = df_tcr_raw[patient_information].copy()
df_patient["Days_from_onset"] = df_patient["Days_from_onset"].astype(
    str
)  # mixed type (str, int)
df_patient = df_patient.drop_duplicates().reset_index(drop=True)

# Assigning barcode as index
df_patient.index = df_patient.pop("barcode")
df_patient.index.name = None
df_patient.head()

adata_tcr.obs[df_patient.columns] = df_patient
adata_tcr.obs.head()

adata_bcr = ir.io.read_10x_vdj(path_bcr_input)
adata_bcr.obs.head(5)

patient_information = ["barcode", "patient_id"]
df_patient = df_bcr_raw[patient_information].copy()
df_patient = df_patient.drop_duplicates().reset_index(drop=True)

# Assigning barcode as index
df_patient.index = df_patient.pop("barcode")
df_patient.index.name = None
df_patient.head()

adata_bcr.obs[df_patient.columns] = df_patient
adata_bcr.obs.head()

ir.tl.chain_qc(adata_tcr)
ir.tl.chain_qc(adata_bcr)

# Visualization
_ = ir.pl.group_abundance(adata_tcr, groupby="Centre", target_col="chain_pairing")
_ = ir.pl.group_abundance(
    adata_tcr, groupby="Centre", target_col="chain_pairing", normalize=True
)
_ = ir.pl.group_abundance(
    adata_tcr, groupby="Centre", target_col="receptor_type", normalize=True
)
_ = ir.pl.group_abundance(
    adata_tcr, groupby="Centre", target_col="receptor_subtype", normalize=True
)
adata_tcr.obs["chain_pairing"].value_counts()

adata_bcr_tmp = adata_bcr[
    adata_bcr.obs["patient_id"].isin(["COVID-003", "AP11"])
].copy()
_ = ir.pl.group_abundance(
    adata_bcr_tmp, groupby="patient_id", target_col="chain_pairing"
)
_ = ir.pl.group_abundance(
    adata_bcr_tmp, groupby="patient_id", target_col="chain_pairing", normalize=True
)
_ = ir.pl.group_abundance(
    adata_bcr_tmp, groupby="patient_id", target_col="receptor_type", normalize=True
)
_ = ir.pl.group_abundance(
    adata_bcr_tmp, groupby="patient_id", target_col="receptor_subtype", normalize=True
)
adata_bcr_tmp.obs["chain_pairing"].value_counts()
print(f"Amount of all B cells:\t\t\t\t{len(adata_bcr)}")
adata_bcr_tmp = adata_bcr[adata_bcr.obs["chain_pairing"] != "no IR"]
print(f"Amount of B cells with AIR:\t\t\t{len(adata_bcr_tmp)}")

adata_bcr_tmp = adata_bcr_tmp[adata_bcr_tmp.obs["chain_pairing"] != "multi_chain"]
print(f"Amount of B cells without dublets:\t\t{len(adata_bcr_tmp)}")

adata_bcr_tmp = adata_bcr_tmp[
    ~adata_bcr_tmp.obs["chain_pairing"].isin(
        ["two full chains", "extra VJ", "extra VDJ"]
    )
]
print(f"Amount of B cells with unique AIR per cell:\t{len(adata_bcr_tmp)}")

adata_bcr_tmp = adata_bcr_tmp[adata_bcr_tmp.obs["chain_pairing"] == "single pair"]
print(f"Amount of B cells with sinlge complete AIR:\t{len(adata_bcr_tmp)}")

adata_bcr_tmp.obs["chain_pairing"].value_counts()
adata_tcr.obs["chain_pairing"].value_counts()

print(f"Amount of all T cells:\t\t\t\t{len(adata_tcr)}")
adata_tcr_tmp = adata_tcr[adata_tcr.obs["chain_pairing"] != "no IR"]
print(f"Amount of T cells with AIR:\t\t\t{len(adata_tcr_tmp)}")

adata_tcr_tmp = adata_tcr_tmp[adata_tcr_tmp.obs["chain_pairing"] != "multi_chain"]
print(f"Amount of T cells without dublets:\t\t{len(adata_tcr_tmp)}")

adata_tcr_tmp = adata_tcr_tmp[
    ~adata_tcr_tmp.obs["chain_pairing"].isin(
        ["two full chains", "extra VJ", "extra VDJ"]
    )
]
print(f"Amount of T cells with unique AIR per cell:\t{len(adata_tcr_tmp)}")

adata_tcr_tmp = adata_tcr_tmp[adata_tcr_tmp.obs["chain_pairing"] == "single pair"]
print(f"Amount of T cells with sinlge complete AIR:\t{len(adata_tcr_tmp)}")

adata_tcr_tmp.obs["chain_pairing"].value_counts()