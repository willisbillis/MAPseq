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
ggsave(paste0("ridgeplot_called_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT, width = OUTPUT_FIG_WIDTH * floor(ncol*0.5))

p = RidgePlot(sc_total,
              features = rownames(sc_total[["HTO"]]),
              ncol = ncol,
              group.by = "HTO_maxID")
ggsave(paste0("ridgeplot_max_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT, width = OUTPUT_FIG_WIDTH * floor(ncol*0.5))

p = RidgePlot(sc_total,
              features = rownames(sc_total[["HTO"]]),
              ncol = ncol,
              group.by = "HTO_classification.global")
ggsave(paste0("ridgeplot_classification_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT, width = OUTPUT_FIG_WIDTH * floor(ncol*0.5))

top_confusion_matrix = as.data.frame(table(sc_total$HTO_maxID[sc_total$HTO_margin > quantile(sc_total$HTO_margin, 0.95)],
       sc_total$HTO_secondID[sc_total$HTO_margin > quantile(sc_total$HTO_margin, 0.95)]))
bot_confusion_matrix = as.data.frame(table(sc_total$HTO_maxID[sc_total$HTO_margin < quantile(sc_total$HTO_margin, 0.05)],
       sc_total$HTO_secondID[sc_total$HTO_margin < quantile(sc_total$HTO_margin, 0.05)]))

top_confusion_matrix$Var1 = as.character(top_confusion_matrix$Var1)
top_confusion_matrix$Var2 = as.character(top_confusion_matrix$Var2)
top_confusion_matrix = top_confusion_matrix[!(top_confusion_matrix$Var1 == top_confusion_matrix$Var2), ]
top_confusion_matrix = top_confusion_matrix[top_confusion_matrix$Freq != 0, ]

bot_confusion_matrix$Var1 = as.character(bot_confusion_matrix$Var1)
bot_confusion_matrix$Var2 = as.character(bot_confusion_matrix$Var2)
bot_confusion_matrix = bot_confusion_matrix[!(bot_confusion_matrix$Var1 == bot_confusion_matrix$Var2), ]
bot_confusion_matrix = bot_confusion_matrix[bot_confusion_matrix$Freq != 0, ]

top_confusion_matrix = top_confusion_matrix[order(top_confusion_matrix$Freq),]
bot_confusion_matrix = bot_confusion_matrix[order(-bot_confusion_matrix$Freq),]
best_htos = c(top_confusion_matrix$Var1[1], top_confusion_matrix$Var2[1])
worst_htos = c(bot_confusion_matrix$Var1[1], bot_confusion_matrix$Var2[1])

Idents(sc_total) = "hash.ID"

p = FeatureScatter(sc_total, cells = colnames(sc_total)[sc_total$hash.ID %in% c(best_htos, "Doublet")],
    feature1 = best_htos[1], feature2 = best_htos[2])
ggsave(paste0("scatter_best.hto.separation_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT, width = OUTPUT_FIG_WIDTH)

p = FeatureScatter(sc_total, cells = colnames(sc_total)[sc_total$hash.ID %in% c(worst_htos, "Doublet")],
    feature1 = worst_htos[1], feature2 = worst_htos[2])
ggsave(paste0("scatter_worst.hto.separation_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT, width = OUTPUT_FIG_WIDTH)

sc_total = subset(sc_total, HTO_classification.global == "Singlet")

DefaultAssay(sc_total) = "RNA"
sc_total <- NormalizeData(sc_total)
sc_total[["percent.Ig"]] <- PercentageFeatureSet(sc_total, pattern = igs)
sc_total[["percent.mt"]] <- PercentageFeatureSet(sc_total, pattern = mt_pattern)

p = DensityScatter(sc_total, "percent.Ig", "percent.mt", quantiles=TRUE)
ggsave("scatter_pct.Ig.v.pct.mt_alldata.png",
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)

p = DensityScatter(sc_total, "nCount_RNA", "nFeature_RNA",
                   quantiles = TRUE, log_x = TRUE, log_y = TRUE)
ggsave("scatter_nCountRNA.v.nFeatRNA_alldata.png",
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
ggsave("scatter_nFeatRNA.v.MT_alldata.png", nFeature_mt_plot,
       width=OUTPUT_FIG_WIDTH, height=OUTPUT_FIG_HEIGHT)
###############################################################################
###############################################################################
# Filtering of cells based on QC criteria

# Expression Filter
stats = data.frame(table(sc_total$patient_id))
colnames(stats) = c("Patients", "Unfiltered_Cells")
rownames(stats) = stats$Patients
stats$Unfiltered_Average_Expression = round(AggregateExpression(sc_total, group.by="patient_id")$RNA %>% colSums / stats$Unfiltered_Cells)

stats$Filtered_by_MT_NF = table(sc_total$patient_id, sc_total$nFeature_RNA <= MIN_GENE_READS & sc_total$percent.mt >= MAX_PCT_MT)[, "TRUE"]
stats$Filtered_by_MT = table(sc_total$patient_id, sc_total$nFeature_RNA <= MIN_GENE_READS)[, "TRUE"] - stats$Filtered_by_MT_NF
stats$Filtered_by_NF = table(sc_total$patient_id, sc_total$percent.mt >= MAX_PCT_MT)[, "TRUE"] - stats$Filtered_by_MT_NF

sc = subset(sc_total,
            subset=percent.mt   < MAX_PCT_MT &
                   nFeature_RNA > MIN_GENE_READS &
                   nFeature_RNA < MAX_GENE_READS)

# HTO Filter
sc <- subset(sc, subset = HTO_classification.global == "Singlet")
sc$hash.ID = factor(sc$hash.ID, levels = unique(sc$hash.ID)) # Remove Doublet and Negative factors

stats$Filtered_Cells = as.numeric(table(sc$patient_id))
stats$Filtered_Average_Expression = round(AggregateExpression(sc, group.by="patient_id")$RNA %>% colSums / stats$Filtered_Cells)

write.csv(stats, "table_patientID_filtering.stats.csv", quote = F, row.names = F)
###############################################################################
# Clustering and DEG (differentially expressed genes) for ADT assay
DefaultAssay(sc) = "RNA"
sc = SCTransform(sc)

# Filter out Ig genes from VariableFeatures, they will clog the results as they are highly-variable by nature
VariableFeatures(sc) = VariableFeatures(sc)[!grepl(igs, VariableFeatures(sc))]

sc = RunPCA(sc, features = VariableFeatures(sc), reduction.name = "rna_pca", reduction.key = "rnapca_", verbose = FALSE)
sc = RunUMAP(sc, reduction = 'rna_pca', dims = 1:40, reduction.key = "rnaumap_", verbose = FALSE)

sc <- FindNeighbors(sc, dims = 1:40, reduction = 'rna_pca', verbose = FALSE)
for (res in c(1, 0.5, 0.25, 0.1, 0.05)) {
  sc = FindClusters(sc, resolution = res, algorithm = 3, verbose = FALSE)
}

p = clustree(sc, prefix="SCT_snn_res.")
ggsave("clustree_clusters_rna.png", p,
       width=OUTPUT_FIG_WIDTH, height=OUTPUT_FIG_HEIGHT)

sc$seurat_clusters_rna = sc[[paste0("SCT_snn_res.", 0.25)]]
Idents(sc) = "seurat_clusters_rna"
sc$seurat_clusters = NULL
sc$seurat_clusters_rna = factor(sc$seurat_clusters_rna)

DefaultAssay(sc) = "RNA"
all_markers = FindAllMarkers(sc)
all_markers = all_markers[all_markers$p_val_adj < 0.05,]
fwrite(all_markers, "DEG_clusters.res0.25.csv", row.names = F, quote = F)

rna_umap = DimPlot(sc, reduction = "umap", group.by = "seurat_clusters_rna")
###############################################################################
# Clustering and DEP (differentially expressed proteins) for ADT assay
DefaultAssay(sc) = "ADT"

sc = NormalizeData(sc)
sc <- ScaleData(sc, features = rownames(sc), verbose = FALSE)
sc <- RunPCA(sc, features = rownames(sc), approx = FALSE, reduction.name = "adt_pca", reduction.key = "adtpca_", verbose = FALSE)
sc = RunUMAP(sc, reduction = 'adt_pca', dims = 1:8, reduction.key = "adtumap_", verbose = FALSE)

sc <- FindNeighbors(sc, dims = 1:8, reduction = 'adt_pca', verbose = FALSE)
for (res in c(1, 0.5, 0.25, 0.1, 0.05)) {
  sc = FindClusters(sc, resolution = res, algorithm = 3, verbose = FALSE)
}

p = clustree(sc, prefix="ADT_snn_res.")
ggsave("clustree_clusters_adt.png", p,
       width=OUTPUT_FIG_WIDTH, height=OUTPUT_FIG_HEIGHT)

sc$seurat_clusters_adt = sc[[paste0("ADT_snn_res.", 0.25)]]
Idents(sc) = "seurat_clusters_adt"
sc$seurat_clusters = NULL
sc$seurat_clusters_adt = factor(sc$seurat_clusters_adt)

all_markers = FindAllMarkers(sc)
all_markers = all_markers[all_markers$p_val_adj < 0.05,]
fwrite(all_markers, "DEP_clusters.res0.25.csv", row.names = F, quote = F)

adt_umap = DimPlot(sc, reduction = "umap", group.by = "seurat_clusters_adt")
###############################################################################
# Creating RNA and ADT UMAP figures
p = rna_umap + adt_umap
ggsave("umap_rna.adt_clusters.png", p, width = OUTPUT_FIG_WIDTH * 2, height = OUTPUT_FIG_HEIGHT)

# save Seurat object
saveRDS(sc, paste0(PROJECT_DIR,"/data/qc_rna.hto.adt_", PROJECT_NAME, ".RDS"))
###############################################################################
###############################################################################
# create ShinyCell app with data - MUST pre-authenticate using shinyapps.io
#      token with rsconnect
if (FALSE) {
  sc = AddMetaData(sc, t(LayerData(sc, assay="HTO")))
  sc = AddMetaData(sc, t(LayerData(sc, assay="ADT")))
  sc_conf = createConfig(sc)
  makeShinyApp(sc, sc_conf, gene.mapping = TRUE,
               shiny.title = paste0(PROJECT_NAME, " RNA + HTO + ADT"),
               shiny.dir = paste0("shiny_",PROJECT_NAME,"_rna"),
               gex.assay="RNA")
}