# run_phase2.QC_LRA_rna.R
# created by M Elliott Williams (https://github.com/willisbillis) Feb 2024

# Install required packages using the package manager 'pacman'
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)
p_load(Seurat, Signac, ggplot2, clustree, dplyr, future, parallel, reticulate,
       multtest, metap, tidyr, scDblFinder, BiocParallel)
p_load_gh("SGDDNB/ShinyCell")
p_load_gh("cellgeni/sceasy")

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("satijalab/azimuth", ref = "master")

# Set python path to ensure reticulate packages can be used
python_path = system("which python", intern = TRUE)
use_python(python_path)
# make sure Seurat v5 is used
options(Seurat.object.assay.version = "v5")
# silence random number warning
options(future.rng.onMisuse = "ignore")
# set seed for reproducibility
set.seed(1234)

## Library descriptions ##
# Seurat: functions for single cell data
# Signac: DensityScatter function for QC
# ggplot2: functions for plotting
# clustree: plotting clusters vs resolution
# dplyr: pipe command '%>%'
# future: multiprocessing limits
# parallel: multiprocessing limits
# reticulate: set which python to use
# multtest: for conserved markers
# metap: for conserved markers
# tidyr: separate function
# Azimuth: cell type annotation
# scDblFinder: doublet detection
# BiocParallel: parallelization for scDblFinder
# ShinyCell: Interact with your data
# sceasy: convert seurat object to anndata format

###############################################################################
#### SET RESOURCE LIMITS ####
###############################################################################
max_cores = 32
max_mem = 128
if (max_cores == -1) {
  max_cores = detectCores()
}
if (max_mem != -1) {
  options("future.globals.maxSize" = (max_mem / max_cores) * 1024^3)
}

plan("multicore", workers = max_cores)
###############################################################################
#### OPTIONS ####
###############################################################################
# REPLACE, must be the same as used in MAPseq pipeline
PROJECT_NAME = "LRA_all_002"
# REPLACE, path to RNA.FB.VDJ analysis dir from MAPseq pipeline
PROJECT_DIR = paste0("/home/Projects/Scharer_sc",
                     "/LRA.MAPseq/LRA_all_002/analysis/RNA.FB.VDJ")
RAW_SEURAT_PATH = paste0(PROJECT_DIR, "/data/raw_rna.hto.adt_",
                         PROJECT_NAME, ".RDS")

GENOME = "hg38"                     # REPLACE, hg38 or mm10
OUTPUT_FIG_WIDTH =  8               # inches, width of output figures
OUTPUT_FIG_HEIGHT = 8               # inches, height of output figures
###############################################################################
#### LOAD DATA AND SET VARIABLES ####
###############################################################################
# TODO: Read a file of Ig instead of pattern matching
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
###############################################################################
#### PLOT DEMULTIPLEXING RESULTS ####
###############################################################################
DefaultAssay(sc_total) = "HTO"
ncol = ceiling(nrow(sc_total[["HTO"]]) / 3)

p = VlnPlot(sc_total,
            features = rownames(sc_total[["HTO"]]),
            ncol = ncol,
            group.by = "HTO_maxID",
            pt.size = 0)
ggsave(paste0("vln_called_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT,
       width = OUTPUT_FIG_WIDTH * floor(ncol * 0.5))

p = VlnPlot(sc_total,
            features = rownames(sc_total[["HTO"]]),
            ncol = ncol,
            group.by = "HTO_classification",
            pt.size = 0)
ggsave(paste0("vln_classification_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT,
       width = OUTPUT_FIG_WIDTH * floor(ncol * 0.5))
###############################################################################
#### CALCULATE QC METRICS (RNA) ####
###############################################################################
# Doublet Detection
DefaultAssay(sc_total) = "RNA"
sc_v3 = sc_total
sc_v3[["HTO"]] = NULL
sc_v3[["ADT"]] = NULL
sc_v3[["RNA"]] = as(sc_v3[["RNA"]], Class = "Assay")
sce <- scDblFinder(as.SingleCellExperiment(sc_v3),
                   samples = "library_id", BPPARAM = MulticoreParam(max_cores))
sc_total$scDblFinder.score <- sce$scDblFinder.score

# Ig and mitochondrial reads detection
sc_total[["percent.Ig"]] <- PercentageFeatureSet(sc_total, pattern = igs)
sc_total[["percent.mt"]] <- PercentageFeatureSet(sc_total, pattern = mt_pattern)

p = DensityScatter(sc_total, "percent.Ig", "percent.mt", quantiles = TRUE)
ggsave("scatter_pct.Ig.v.pct.mt_alldata.png",
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)

p = DensityScatter(sc_total, "scDblFinder.score", "percent.mt",
                   quantiles = TRUE)
ggsave("scatter_dbl.v.pct.mt_alldata.png",
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)

p = DensityScatter(sc_total, "nCount_RNA", "nFeature_RNA",
                   quantiles = TRUE, log_x = TRUE, log_y = TRUE)
ggsave("scatter_nCountRNA.v.nFeatRNA_alldata.png",
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)

p = DensityScatter(sc_total, "nCount_HTO", "nCount_RNA",
                   quantiles = TRUE, log_x = TRUE, log_y = TRUE)
ggsave("scatter_nCountHTO.v.nCountRNA_alldata.png",
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)

p = DensityScatter(sc_total, "nFeature_HTO", "nFeature_RNA",
                   quantiles = TRUE, log_x = TRUE, log_y = TRUE)
ggsave("scatter_nFeatHTO.v.nFeatRNA_alldata.png",
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)
###############################################################################
#### RNA CUTOFFS ####
###############################################################################
# PAUSE, view scatter figures above and determine appropriate cutoffs below
MAX_PCT_MT = 5        # REPLACE, maximum percent mitochondrial reads per cell
DBL_LIMIT = 0.8       # REPLACE, maximum scDblFinder score to permit
MIN_GENE_READS = 100   # REPLACE, minimum genes with reads per cell
MAX_GENE_READS = Inf  # REPLACE, maximum genes with reads per cell
#                                (set plasma cell limit to Inf)

p = DensityScatter(sc_total, "nFeature_RNA", "percent.mt",
                   quantiles=TRUE, log_x=TRUE)
nfeat_mt_plot <- p +
  geom_hline(yintercept = MAX_PCT_MT, linetype = "dashed") +
  geom_vline(xintercept = c(MIN_GENE_READS, MAX_GENE_READS),
             linetype = "dashed")
ggsave("scatter_nFeatRNA.v.pct.mt_filtered.png", nfeat_mt_plot,
       width=OUTPUT_FIG_WIDTH, height=OUTPUT_FIG_HEIGHT)
###############################################################################
#### QUANTIFY QC FILTERING ####
###############################################################################
# adjust metadata to accomodate Seurat's AggregateExpression
sc_total$library_id = gsub("_", "-", sc_total$library_id)
sc_total$patient_id = gsub("_", "-", sc_total$patient_id)
sc_total$patient_id[sc_total$HTO_maxID == "Doublet"] = "Doublet"
sc_total$patient_id[sc_total$HTO_maxID == "Negative"] = "Negative"
sc_total = sc_total[, !is.na(sc_total$patient_id)]
# create new column for unique sample ID - adjust as needed for each dataset
sc_total$sample_id = paste(sc_total$library_id,
                           sc_total$patient_id,
                           sep = "-")
stats = data.frame(sample_id = unique(sc_total$sample_id))
sample_id_counts = as.data.frame(table(sc_total$sample_id))
stats$Unfiltered_Cells = sample_id_counts$Freq[match(stats$sample_id,
                                                     sample_id_counts$Var1)]
stats[is.na(stats)] = 0
cells_rna = AggregateExpression(sc_total,
                                group.by = "sample_id")$RNA %>%
  colSums
cells_adt = AggregateExpression(sc_total,
                                group.by = "sample_id")$ADT %>%
  colSums
cells_hto = AggregateExpression(sc_total,
                                group.by = "sample_id")$HTO %>%
  colSums

stats$Unfiltered_Avg_Expression.RNA = round(cells_rna[match(stats$sample_id,
                                                            names(cells_rna))] /
                                              stats$Unfiltered_Cells)
stats$Unfiltered_Avg_Expression.ADT = round(cells_adt[match(stats$sample_id,
                                                            names(cells_adt))] /
                                              stats$Unfiltered_Cells)
stats$Unfiltered_Avg_Expression.HTO = round(cells_hto[match(stats$sample_id,
                                                            names(cells_hto))] /
                                              stats$Unfiltered_Cells)

sc = subset(sc_total,
            subset = percent.mt < MAX_PCT_MT &
              scDblFinder.score < DBL_LIMIT &
              nFeature_RNA > MIN_GENE_READS &
              nFeature_RNA < MAX_GENE_READS &
              HTO_maxID != "Doublet")

sample_id_counts = as.data.frame(table(sc$sample_id))
stats$Filtered_Cells = sample_id_counts$Freq[match(stats$sample_id,
                                                   sample_id_counts$Var1)]
stats[is.na(stats)] = 0

cells_rna = AggregateExpression(sc, group.by = "sample_id")$RNA %>% colSums
cells_adt = AggregateExpression(sc, group.by = "sample_id")$ADT %>% colSums
cells_hto = AggregateExpression(sc, group.by = "sample_id")$HTO %>% colSums

stats$Filtered_Avg_Expression.RNA = round(cells_rna[match(stats$sample_id,
                                                          names(cells_rna))] /
                                            stats$Filtered_Cells)
stats$Filtered_Avg_Expression.ADT = round(cells_adt[match(stats$sample_id,
                                                          names(cells_adt))] /
                                            stats$Filtered_Cells)
stats$Filtered_Avg_Expression.HTO = round(cells_hto[match(stats$sample_id,
                                                          names(cells_hto))] /
                                            stats$Filtered_Cells)
stats[is.na(stats)] = 0
write.csv(stats, "QC.rna_sampleID_filtering.stats.csv",
          quote = FALSE, row.names = FALSE)
###############################################################################
# SAVE RAW SEURAT OBJECT
###############################################################################
saveRDS(sc_total,
        paste0(PROJECT_DIR, "/data/raw_rna.hto.adt_", PROJECT_NAME, ".RDS"))
###############################################################################
#### BATCH CORRECTION (OPTIONAL) ####
###############################################################################
if (FALSE) {
  #### RNA Harmony Batch correction ####
  DefaultAssay(sc) = "RNA"
  split_column = "sample_id"
  sc[["RNA"]] = split(sc[["RNA"]], f = sc@meta.data[[split_column]], layers = "counts")
  sc <- SCTransform(sc, verbose = FALSE)

  # Filter out Ig genes from VariableFeatures, they will clog the results as
  #     they are highly-variable by nature
  non_ig_mask = !grepl(igs, VariableFeatures(sc))
  VariableFeatures(sc) = VariableFeatures(sc)[non_ig_mask]
  sc <- RunPCA(sc, npcs = 15,
               reduction.name = "rna.pca", verbose = FALSE)
  sc <- FindNeighbors(sc, dims = 1:15, reduction = "rna.pca",
                      verbose = FALSE)
  sc <- FindClusters(sc, resolution = 2, algorithm = 4,
                     cluster.name = "unintegrated_rna.clusters",
                     verbose = FALSE)
  sc <- RunUMAP(sc, dims = 1:15, reduction = "rna.pca",
                reduction.name = "umap.rna.unintegrated",
                verbose = FALSE)
  # visualize by batch annotations
  p = DimPlot(sc, reduction = "umap.rna.unintegrated", group.by = split_column)
  ggsave(paste0("umap_rna.unintegrated_", split_column, ".pdf"), p, width = 8, height = 6)
  ggsave(paste0("umap_rna.unintegrated_", split_column, ".png"), p, width = 8, height = 6)

  sc <- IntegrateLayers(
    object = sc, method = HarmonyIntegration,
    orig = "rna.pca", new.reduction = "integrated.rna.harmony",
    verbose = FALSE
  )
  sc = JoinLayers(sc, assay = "RNA")

  sc <- FindNeighbors(sc, reduction = "integrated.rna.harmony",
                      dims = 1:15, verbose = FALSE)
  sc <- FindClusters(sc, resolution = 2, algorithm = 4,
                     cluster.name = "rna.clusters",
                     verbose = FALSE)
  sc <- RunUMAP(sc, reduction = "integrated.rna.harmony",
                dims = 1:15, reduction.name = "umap.rna",
                verbose = FALSE)
  p <- DimPlot(sc, reduction = "integrated.rna.harmony",
               group.by = split_column)
  ggsave(paste0("umap_rna.integrated_", split_column, ".pdf"), p, width = 8, height = 6)
  ggsave(paste0("umap_rna.integrated_", split_column, ".png"), p, width = 8, height = 6)

  #### ADT Harmony Batch correction ####
  DefaultAssay(sc) = "ADT"
  sc[["ADT"]] = split(sc[["ADT"]], f = sc[[split_column]])
  VariableFeatures(sc) <- rownames(sc[["ADT"]])
  sc <- NormalizeData(sc, normalization.method = "CLR", margin = 2,
                      verbose = FALSE)
  sc <- ScaleData(sc, features = rownames(sc[["ADT"]]),
                  do.center = TRUE,
                  do.scale = FALSE,
                  verbose = FALSE)
  sc <- RunPCA(sc, npcs = 15,
               reduction.name = "adt.pca", verbose = FALSE)
  sc <- FindNeighbors(sc, dims = 1:15, reduction = "adt.pca",
                      verbose = FALSE)
  sc <- FindClusters(sc, resolution = 2, algorithm = 4,
                     cluster.name = "unintegrated_adt.clusters",
                     verbose = FALSE)
  sc <- RunUMAP(sc, dims = 1:15, reduction = "adt.pca",
                reduction.name = "umap.adt.unintegrated",
                verbose = FALSE)
  # visualize by batch annotations
  p = DimPlot(sc, reduction = "umap.adt.unintegrated", group.by = split_column)
  ggsave(paste0("umap_adt.unintegrated_", split_column, ".pdf"), p, width = 8, height = 6)
  ggsave(paste0("umap_adt.unintegrated_", split_column, ".png"), p, width = 8, height = 6)

  sc <- IntegrateLayers(
    object = sc, method = HarmonyIntegration,
    features = rownames(sc[["ADT"]]),
    orig = "adt.pca", new.reduction = "integrated.adt.harmony",
    verbose = FALSE
  )
  sc = JoinLayers(sc)

  sc <- FindNeighbors(sc, reduction = "integrated.adt.harmony",
                      dims = 1:15, verbose = FALSE)
  sc <- FindClusters(sc, resolution = 2, algorithm = 4,
                     cluster.name = "adt.clusters",
                     verbose = FALSE)
  sc <- RunUMAP(sc, reduction = "integrated.adt.harmony",
                dims = 1:15, reduction.name = "umap.adt",
                verbose = FALSE)
  p <- DimPlot(sc, reduction = "integrated.adt.harmony",
               group.by = split_column)
  ggsave(paste0("umap_adt.integrated_", split_column, ".pdf"), p, width = 8, height = 6)
  ggsave(paste0("umap_adt.integrated_", split_column, ".png"), p, width = 8, height = 6)

  #### (Optional) WNN Integration of RNA and ADT assays ####
  sc <- FindMultiModalNeighbors(
    sc, reduction.list = list("rna.pca", "adt.pca"),
    dims.list = list(1:15, 1:15), verbose = FALSE
  )
  sc = RunUMAP(sc, nn.name = "weighted.nn",  reduction.name = "wnn.umap",
               reduction.key = "wnnUMAP_", verbose = FALSE)
  p <- DimPlot(sc, reduction = "wnn.umap", group.by = split_column)
  ggsave(paste0("umap_rna.adt.wnn_", split_column, ".pdf"), p, width = 8, height = 6)
  ggsave(paste0("umap_rna.adt.wnn_", split_column, ".png"), p, width = 8, height = 6)
}
###############################################################################
#### NON-BATCH CORRECTED DIMENSIONAL REDUCTION ####
###############################################################################
DefaultAssay(sc) = "RNA"
sc <- SCTransform(sc, verbose = FALSE)
# Filter out Ig genes from VariableFeatures, they will clog the results as
#     they are highly-variable by nature
non_ig_mask = !grepl(igs, VariableFeatures(sc))
VariableFeatures(sc) = VariableFeatures(sc)[non_ig_mask]
sc <- RunPCA(sc, npcs = 40,
             reduction.name = "rna.pca", verbose = FALSE)
sc <- FindNeighbors(sc, dims = 1:40, reduction = "rna.pca",
                    verbose = FALSE)
sc <- RunUMAP(sc, dims = 1:40, reduction = "rna.pca",
              reduction.name = "umap.rna",
              verbose = FALSE)

DefaultAssay(sc) = "ADT"
VariableFeatures(sc) <- rownames(sc[["ADT"]])
sc = NormalizeData(sc, normalization.method = "CLR", margin = 2)
sc = ScaleData(sc)
###############################################################################
#### CLUSTERING AND ANNOTATION ####
###############################################################################
# Annotate PBMC cell types using Azimuth's PBMC reference
DefaultAssay(sc) = "SCT"
sc_v3 = sc
sc_v3[["ADT"]] = NULL
sc_v3[["RNA"]] = NULL
sc_v3[["SCT"]] = as(sc_v3[["SCT"]], Class = "Assay")
# REPLACE AZIMUTH REFERENCE WITH APPROPRIATE DATASET
sc_v3 <- RunAzimuth(sc_v3, reference = "pbmcref")

sc$predicted.celltype.l1 = sc_v3$predicted.celltype.l1
sc$predicted.celltype.l2 = sc_v3$predicted.celltype.l2
sc$predicted.celltype.l3 = sc_v3$predicted.celltype.l3
sc$predicted.celltype.l1.score = sc_v3$predicted.celltype.l1.score
sc$predicted.celltype.l2.score = sc_v3$predicted.celltype.l2.score
sc$predicted.celltype.l3.score = sc_v3$predicted.celltype.l3.score

sc@reductions$ref.umap = sc_v3@reductions$ref.umap

# Different cluster resolutions for SCT
graph = "SCT_snn"
for (res in c(1, 0.5, 0.25, 0.1, 0.05)) {
  sc = FindClusters(sc, resolution = res, graph.name = graph,
                    algorithm = 3, verbose = FALSE)
}

p = clustree(sc, prefix = paste0(graph, "_res."))
ggsave("clustree_clusters_rna.png", p,
       width=OUTPUT_FIG_WIDTH, height=OUTPUT_FIG_HEIGHT)

sc$seurat_clusters = sc[[paste0(graph, "_res.", 0.25)]]
Idents(sc) = "seurat_clusters"
sc$seurat_clusters = factor(sc$seurat_clusters)

# DEG testing between clusters (one vs all)
all_markers = FindAllMarkers(sc, verbose = FALSE, assay = "SCT")
all_markers = all_markers[all_markers$p_val_adj < 0.05, ]
write.csv(all_markers, paste("DEG_", graph, ".clusters.res0.25.csv"),
          row.names = FALSE, quote = FALSE)

# DEP testing between clusters (one vs all)
all_markers = FindAllMarkers(sc, verbose = FALSE, assay = "ADT")
all_markers = all_markers[all_markers$p_val_adj < 0.05, ]
write.csv(all_markers, paste("DEP_", graph, ".clusters.res0.25.csv"),
          row.names = FALSE, quote = FALSE)
###############################################################################
# SAVE SEURAT OBJECT AND SESSION INFO
###############################################################################
saveRDS(sc, paste0(PROJECT_DIR, "/data/qc_rna.hto.adt_", PROJECT_NAME, ".RDS"))
# Save the R session environment information
capture.output(sessionInfo(),
               file = paste0(PROJECT_DIR, "/",
                             PROJECT_NAME,
                             ".Rsession.Info.",
                             gsub("\\D", "", Sys.time()), ".txt"))
###############################################################################
# Optional Additonal Analyses
###############################################################################
# create ShinyCell app with data - MUST pre-authenticate using shinyapps.io
#      token with rsconnect
if (FALSE) {
  adt_cts = LayerData(sc, assay = "ADT", layer = "counts")
  sct_cts = LayerData(sc, assay = "SCT", layer = "counts")
  adt_data = LayerData(sc, assay = "ADT", layer = "data")
  sct_data = LayerData(sc, assay = "SCT", layer = "data")
  sct_act_cts = rbind(sct_cts, adt_cts)
  sct_act_data = rbind(sct_data, adt_data)
  sc[["SCT_ADT"]] = CreateAssay5Object(counts = sct_act_cts,
                                       data = sct_act_data)
  DefaultAssay(sc) = "SCT_ADT"
  sc = FindVariableFeatures(sc)
  sc_conf = createConfig(sc)
  makeShinyApp(sc, sc_conf, gene.mapping = FALSE,
               shiny.title = paste0(PROJECT_NAME, " RNA + ADT + HTO"),
               shiny.dir = paste0("shiny_", PROJECT_NAME, "_rna"),
               gex.assay = "SCT_ADT")
  rsconnect::deployApp(paste0("shiny_", PROJECT_NAME, "_rna"))
}
###############################################################################
# Save h5ad for CellxGene use (https://github.com/chanzuckerberg/cellxgene)
if (FALSE) {
  adt_cts = LayerData(sc, assay="ADT", layer = "counts")
  sct_cts = LayerData(sc, assay="SCT", layer = "counts")
  adt_data = LayerData(sc, assay="ADT", layer = "data")
  sct_data = LayerData(sc, assay="SCT", layer = "data")
  sct_act_cts = rbind(sct_cts, adt_cts)
  sct_act_data = rbind(sct_data, adt_data)
  sc[["SCT_ADT"]] = CreateAssay5Object(counts = sct_act_cts,
                                       data = sct_act_data)
  DefaultAssay(sc) = "SCT_ADT"
  sc[["SCT_ADT"]] = as(sc[["SCT_ADT"]], Class = "Assay")
  sc[["HTO"]] = as(sc[["HTO"]], Class = "Assay")
  sc[["RNA"]] = as(sc[["RNA"]], Class = "Assay")
  sc[["SCT"]] = as(sc[["SCT"]], Class = "Assay")
  sc[["HTO"]] = as(sc[["HTO"]], Class = "Assay")
  sceasy::convertFormat(sc, from = "seurat", to = "anndata",
                        outFile = paste0("qc_sct.adt_", PROJECT_NAME, ".h5ad"))
}