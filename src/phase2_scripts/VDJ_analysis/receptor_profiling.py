import pandas as pd
import scirpy as ir

# Ignore DtypeWarnings from pandas
warnings.simplefilter(action="ignore", category=pd.errors.DtypeWarning)

path_data = "data/"
path_bcr_input = f"{path_data}/BCR_00_read_aligned.csv"
path_tcr_input = f"{path_data}/TCR_00_read_aligned.tsv"

path_bcr_out = f"{path_data}/BCR_01_preprocessed.h5ad" # This file is never written to.

path_tcr_csv = f"{path_data}/TCR_00_read_aligned.csv" # Temporary CSV conversion
path_tcr_out = f"{path_data}/TCR_01_preprocessed.h5ad" # This file is never written to.


# Load BCR data
df_bcr_raw = pd.read_csv(path_bcr_input, index_col=0)

# Ensure 'productive' column is string type
df_bcr_raw["productive"] = df_bcr_raw["productive"].astype(str)
print(f"Total BCR measurements: {len(df_bcr_raw)}")


# Load TCR data (from TSV)
df_tcr_raw = pd.read_csv(path_tcr_input, sep="\t")
df_tcr_raw["barcode"] = df_tcr_raw.pop("CellID")  # Rename CellID to barcode

# Save TCR data as CSV (for scirpy compatibility)
df_tcr_raw.to_csv(path_tcr_csv, index=False) # Add index=False to prevent writing row numbers
print(f"Total TCR measurements: {len(df_tcr_raw)}")


# Read TCR data using scirpy
adata_tcr = ir.io.read_10x_vdj(path_tcr_csv)
print(f"Number of T cells: {len(adata_tcr)}")


# Extract and clean patient information from TCR data
patient_information = [
    "barcode", "Centre", "Sample", "patient_id", "Collection_Day", "Sex", "Swab_result",
    "Status", "Smoker", "Status_on_day_collection", "Status_on_day_collection_summary",
    "Days_from_onset", "time_after_LPS", "Worst_Clinical_Status", "Outcome",
    "initial_clustering", "study_id", "AgeRange", "Age"
]
df_patient = df_tcr_raw[patient_information].copy()
df_patient["Days_from_onset"] = df_patient["Days_from_onset"].astype(str)
df_patient = df_patient.drop_duplicates().reset_index(drop=True)
df_patient.set_index("barcode", inplace=True) # Set index directly
adata_tcr.obs = adata_tcr.obs.join(df_patient, how="left") # Efficiently join dataframes


# Read BCR data using scirpy
adata_bcr = ir.io.read_10x_vdj(path_bcr_input)
print(f"Number of B cells: {len(adata_bcr)}")

# Extract patient information from BCR data
patient_information = ["barcode", "patient_id"]
df_patient = df_bcr_raw[patient_information].copy()
df_patient = df_patient.drop_duplicates().reset_index(drop=True)
df_patient.set_index("barcode", inplace=True)
adata_bcr.obs = adata_bcr.obs.join(df_patient, how="left")

# Perform chain QC
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

# Save the processed data
adata_tcr.write_h5ad(path_tcr_out)
adata_bcr.write_h5ad(path_bcr_out)