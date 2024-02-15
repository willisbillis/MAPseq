# run_seuratQC_mapseq.R - created by M Elliott Williams (https://github.com/willisbillis) Feb 2024

# Install required packages using the package manager 'pacman'
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)
p_load(Seurat,Signac,GenomeInfoDb,
    EnsDb.Mmusculus.v79,EnsDb.Hsapiens.v86,
    ggplot2,data.table,clustree,tidyr,dplyr,
    tidytext)
p_load_gh("satijalab/seurat-data")
p_load_gh("mojaveazure/seurat-disk")
p_load_gh("SGDDNB/ShinyCell")

set.seed(1234)                # set seed for reproducibility
## Library descriptions ##
# Seurat: functions for single cell data
# Signac: DensityScatter function for QC
# SeuratData: save Seurat object
# SeuratDisk: save Seurat object
# GenomeInfoDb: database for genomic annotations
# EnsDb.Mmusculus.v79: database for mm10 annotations
# EnsDb.Hsapiens.v86: database for hg38 annotations
# ggplot2: functions for plotting
# data.table: write out csv quickly
# clustree: plotting clusters vs resolution
# tidyr: for "separate" function
# dplyr: for "slice_max" function
# tidytext: reorder_within()
# ShinyCell: Interact with your data

###############################################################################
## OPTIONS

# REPLACE, must be the same as used in MAPseq pipeline
PROJECT_NAME = "LRA_all"
# REPLACE, path to output dir from MAPseq pipeline
PROJECT_DIR = "/home/boss_lab/Projects/Scharer_sc/LRA.MAPseq/LRA_all/analysis/RNA.FB.VDJ"
RAW_SEURAT_PATH = paste0(PROJECT_DIR,"/data/raw_rna.hto.adt_", PROJECT_NAME, ".RDS")

GENOME = "hg38"                     # REPLACE, hg38 or mm10
OUTPUT_FIG_WIDTH =  8               # inches, width of output figures
OUTPUT_FIG_HEIGHT = 8

###############################################################################
# TODO: Read a file of Ig instead of pattern matching?
if (GENOME == "hg38") {
  mt_pattern = "^MT-"
  igs = "^IGH|^IGK|^IGL"
} else if (GENOME == "mm10") {
  mt_pattern = "^mt-"
  igs = "^Igh|^Igk|^Igl"
} else {
  print(paste("[ERROR] Requested genome", GENOME,
              "is not available or not implemented. Exiting."))
  quit()
}

setwd(PROJECT_DIR)
sc_total = readRDS(RAW_SEURAT_PATH)

DefaultAssay(sc_total) = "HTO"
ncol = ceiling(nrow(sc_total[["HTO"]]) / 3)

p = RidgePlot(sc_total,
              features = rownames(sc_total[["HTO"]]),
              ncol = ncol,
              group.by = "hash.ID")
ggsave(paste0(PROJECT_DIR, "/ridgeplot_called_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT, width = OUTPUT_FIG_WIDTH * floor(ncol*0.5))

p = RidgePlot(sc_total,
              features = rownames(sc_total[["HTO"]]),
              ncol = ncol,
              group.by = "HTO_maxID")
ggsave(paste0(PROJECT_DIR, "/ridgeplot_max_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT, width = OUTPUT_FIG_WIDTH * floor(ncol*0.5))

p = RidgePlot(sc_total,
              features = rownames(sc_total[["HTO"]]),
              ncol = ncol,
              group.by = "HTO_classification.global")
ggsave(paste0(PROJECT_DIR, "/ridgeplot_classification_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT, width = OUTPUT_FIG_WIDTH * floor(ncol*0.5))

confusion_matrix = as.data.frame(table(sc_total$HTO_maxID, sc_total$HTO_secondID))
confusion_matrix = confusion_matrix[!(confusion_matrix$Var1 == confusion_matrix$Var2), ]
confusion_matrix = confusion_matrix[confusion_matrix$Freq != 0, ]
confusion_matrix$Var1 = as.character(confusion_matrix$Var1)
confusion_matrix$Var2 = as.character(confusion_matrix$Var2)

min_idx = as.numeric(rownames(confusion_matrix)[confusion_matrix$Freq == min(confusion_matrix$Freq)][1])
max_idx = as.numeric(rownames(confusion_matrix)[confusion_matrix$Freq == max(confusion_matrix$Freq)])

Idents(sc_total) = "hash.ID"

p = FeatureScatter(sc_total,
    feature1 = confusion_matrix$Var1[rownames(confusion_matrix) == min_idx], feature2 = confusion_matrix$Var2[rownames(confusion_matrix) == min_idx])
ggsave(paste0(PROJECT_DIR, "/scatter_good.hto.separation_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT, width = OUTPUT_FIG_WIDTH)

p = FeatureScatter(sc_total,
    feature1 = confusion_matrix$Var1[rownames(confusion_matrix) == max_idx], feature2 = confusion_matrix$Var2[rownames(confusion_matrix) == max_idx])
ggsave(paste0(PROJECT_DIR, "/scatter_worst.hto.separation_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT, width = OUTPUT_FIG_WIDTH)

sc_total = subset(sc_total, HTO_classification.global == "Singlet")

DefaultAssay(sc_total) = "RNA"
sc_total <- NormalizeData(sc_total)
sc_total[["percent.Ig"]] <- PercentageFeatureSet(sc_total, pattern = igs)
sc_total[["percent.mt"]] <- PercentageFeatureSet(sc_total, pattern = mt_pattern)

p = DensityScatter(sc_total, "percent.Ig", "percent.mt", quantiles=TRUE)
ggsave(paste0(PROJECT_DIR, "/scatter_pct.Ig.v.pct.mt_alldata.png"),
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)

p = DensityScatter(sc_total, "nCount_RNA", "nFeature_RNA",
                   quantiles = TRUE, log_x = TRUE, log_y = TRUE)
ggsave(paste0(PROJECT_DIR, "/scatter_nCountRNA.v.nFeatRNA_alldata.png"),
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)
###############################################################################
###############################################################################
# PAUSE, view scatter figures above and determine appropriate cutoffs below
MAX_PCT_MT = 5        # REPLACE, maximum percent mitochondrial reads per cell
MIN_GENE_READS = 300   # REPLACE, minimum genes with reads per cell
MAX_GENE_READS = Inf  # REPLACE, maximum genes with reads per cell; set plasma cell limit to Inf

p = DensityScatter(sc_total, "nFeature_RNA", "percent.mt", quantiles=TRUE, log_x=TRUE)
nFeature_mt_plot <- p +
  geom_hline(yintercept=MAX_PCT_MT,linetype='dashed') +
  geom_vline(xintercept=c(MIN_GENE_READS, MAX_GENE_READS), linetype='dashed')
ggsave(paste0(PROJECT_DIR, "/scatter_nFeatRNA.v.MT_alldata.png"), nFeature_mt_plot,
       width=OUTPUT_FIG_WIDTH, height=OUTPUT_FIG_HEIGHT)
###############################################################################
###############################################################################
categorize_cells = function(sc_seurat, cat_column, names_prefix) {
  table(sc_seurat$patient_id, unlist(sc_seurat[[cat_column]]), dnn = list("Patients", cat_column)) %>%
    data.frame %>%
    pivot_wider(names_from = all_of(cat_column), values_from = "Freq", names_prefix = names_prefix) %>%
    select(!Patients)
}

# Expression Filter
stats = data.frame(table(sc_total$patient_id))
colnames(stats) = c("Patients", "Unfiltered_Cells")
rownames(stats) = stats$Patients
stats$Unfiltered_Average_Expression = AggregateExpression(sc_total, group.by="patient_id")$RNA %>% colSums / stats$Unfiltered_Cells
stats = cbind(stats, categorize_cells(sc_total, "hash.ID", "Unfiltered_"))

stats$Filtered_by_MT_NF = table(sc_total$patient_id, sc_total$nFeature_RNA <= MIN_GENE_READS & sc_total$percent.mt >= MAX_PCT_MT)[, "TRUE"]
stats$Filtered_by_MT = table(sc_total$patient_id, sc_total$nFeature_RNA <= MIN_GENE_READS)[, "TRUE"] - stats$Filtered_by_MT_NF
stats$Filtered_by_NF = table(sc_total$patient_id, sc_total$percent.mt >= MAX_PCT_MT)[, "TRUE"] - stats$Filtered_by_MT_NF

sc = subset(sc_total,
            subset=percent.mt   < MAX_PCT_MT &
                   nFeature_RNA > MIN_GENE_READS &
                   nFeature_RNA < MAX_GENE_READS)

stats$Filtered_Cells = as.numeric(table(sc$patient_id))
stats = cbind(stats, categorize_cells(sc_total, "hash.ID", "Filtered_"))

# HTO Filter
sc <- subset(sc, subset = HTO_classification.global == "Singlet")
sc$hash.ID = factor(sc$hash.ID, levels = unique(sc$hash.ID)) # Remove Doublet and Negative factors

stats$Filtered_Singleton_Cells = as.numeric(table(sc$patient_id))
stats$Filtered_Singleton_Average_Expression = AggregateExpression(sc, group.by="patient_id")$RNA %>% colSums %>% round

write.csv(stats, "table_patientID_filtering.stats.csv", quote = F, row.names = F)