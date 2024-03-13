# run_seuratQC_mapseq_rna.R - created by M Elliott Williams (https://github.com/willisbillis) Feb 2024

# Install required packages using the package manager 'pacman'
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)
p_load(Seurat, Signac, ggplot2, clustree, dplyr, future, parallel, reticulate,
       multtest, metap)
p_load_gh("SGDDNB/ShinyCell")
p_load_gh("cellgeni/sceasy")

# Set python path to ensure reticulate packages can be used
python_path = system("which python", intern = TRUE)
use_python(python_path)
# make sure you are running Seurat v5
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
# ShinyCell: Interact with your data

###############################################################################
#### SET RESOURCE LIMITS ####
max_cores = 32
max_mem = 32
if (max_cores == -1) {
  max_cores = detectCores()
}
if (max_mem != -1) {
  options("future.globals.maxSize" = (max_mem / max_cores) * 1024^3)
}

plan("multicore", workers = max_cores)
###############################################################################
## OPTIONS

# REPLACE, must be the same as used in MAPseq pipeline
PROJECT_NAME = "LRA_all"
# REPLACE, path to RNA.FB.VDJ analysis dir from MAPseq pipeline
PROJECT_DIR = "/home/boss_lab/Projects/Scharer_sc/LRA.MAPseq/LRA_all/analysis/RNA.FB.VDJ"
RAW_SEURAT_PATH = paste0(PROJECT_DIR,"/data/raw_rna.hto.adt_", PROJECT_NAME, ".RDS")
HTO_DEMUX_PATH = paste0(PROJECT_DIR,"/../../pipeline/RNA.FB.VDJ/hashtag_ref_rna.csv")

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

p = VlnPlot(sc_total,
            features = rownames(sc_total[["HTO"]]),
            ncol = ncol,
            group.by = "hash.ID",
            pt.size = 0)
ggsave(paste0("vln_called_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT, width = OUTPUT_FIG_WIDTH * floor(ncol*0.5))

p = VlnPlot(sc_total,
            features = rownames(sc_total[["HTO"]]),
            ncol = ncol,
            group.by = "HTO_maxID",
            pt.size = 0)
ggsave(paste0("vln_max_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT, width = OUTPUT_FIG_WIDTH * floor(ncol*0.5))

p = VlnPlot(sc_total,
            features = rownames(sc_total[["HTO"]]),
            ncol = ncol,
            group.by = "HTO_classification.global",
            pt.size = 0)
ggsave(paste0("vln_classification_", PROJECT_NAME, ".png"),
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

confusion_matrix_all = rbind(top_confusion_matrix, bot_confusion_matrix)
colnames(confusion_matrix_all) = c("HT_1st", "HT_2nd", "mixing_degree")
write.csv(confusion_matrix_all, paste0("HTC.combos_",PROJECT_NAME,"_metrics.csv"), quote = F, row.names = F)

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
ggsave("scatter_nFeatRNA.v.MT_filtered.png", nFeature_mt_plot,
       width=OUTPUT_FIG_WIDTH, height=OUTPUT_FIG_HEIGHT)
###############################################################################
###############################################################################
# Filtering of cells based on QC criteria

# Expression Filter
hto_reference = read.csv(HTO_DEMUX_PATH)
stats = data.frame(Patients = hto_reference$patient_id)
patient_id.counts = as.data.frame(table(sc_total$patient_id))
stats$Unfiltered_Cells = patient_id.counts$Freq[match(stats$Patients, patient_id.counts$Var1)]
stats[is.na(stats)] = 0
cell_counts.rna = round(AggregateExpression(sc_total, group.by="patient_id")$RNA %>% colSums)
cell_counts.adt = round(AggregateExpression(sc_total, group.by="patient_id")$ADT %>% colSums)
cell_counts.hto = round(AggregateExpression(sc_total, group.by="patient_id")$HTO %>% colSums)

stats$Unfiltered_Average_Expression.RNA = cell_counts.rna[match(stats$Patients, names(cell_counts.rna))] /
  stats$Unfiltered_Cells
stats$Unfiltered_Average_Expression.ADT = cell_counts.adt[match(stats$Patients, names(cell_counts.adt))] /
  stats$Unfiltered_Cells
stats$Unfiltered_Average_Expression.HTO = cell_counts.hto[match(stats$Patients, names(cell_counts.hto))] /
  stats$Unfiltered_Cells

sc = subset(sc_total,
            subset = percent.mt < MAX_PCT_MT &
              nFeature_RNA > MIN_GENE_READS &
              nFeature_RNA < MAX_GENE_READS)

# HTO Filter
sc <- subset(sc, subset = HTO_classification.global != "Doublet")
DefaultAssay(sc) = "ADT"
VariableFeatures(sc) <- rownames(sc[["ADT"]])
sc = NormalizeData(sc, normalization.method = 'CLR', margin = 2)
sc = ScaleData(sc)

patient_id.counts = as.data.frame(table(sc$patient_id))
stats$Filtered_Cells = patient_id.counts$Freq[match(stats$Patients, patient_id.counts$Var1)]
stats[is.na(stats)] = 0

cell_counts.rna = round(AggregateExpression(sc, group.by="patient_id")$RNA %>% colSums)
cell_counts.adt = round(AggregateExpression(sc, group.by="patient_id")$ADT %>% colSums)
cell_counts.hto = round(AggregateExpression(sc, group.by="patient_id")$HTO %>% colSums)

stats$Filtered_Average_Expression.RNA = cell_counts.rna[match(stats$Patients, names(cell_counts.rna))] /
  stats$Filtered_Cells
stats$Filtered_Average_Expression.ADT = cell_counts.adt[match(stats$Patients, names(cell_counts.adt))] /
  stats$Filtered_Cells
stats$Filtered_Average_Expression.HTO = cell_counts.hto[match(stats$Patients, names(cell_counts.hto))] /
  stats$Filtered_Cells
stats[is.na(stats)] = 0

write.csv(stats, "QC.rna_patientID_filtering.stats.csv", quote = F, row.names = F)
###############################################################################
# Batch correction, Clustering, DEG (differentially expressed genes), and
#     DEP (differentially expressed proteins) for RNA assay

#### RNA RPCA Batch correction ####
DefaultAssay(sc) = "RNA"
sc[["RNA"]] = split(sc[["RNA"]], f = sc$patient_id)
sc <- SCTransform(sc, verbose = FALSE)

# Filter out Ig genes from VariableFeatures, they will clog the results as
#     they are highly-variable by nature
non_ig_mask = !grepl(igs, VariableFeatures(sc))
VariableFeatures(sc) = VariableFeatures(sc)[non_ig_mask]
sc <- RunPCA(sc, npcs = 40,
             reduction.name = "rna.pca", verbose = FALSE)
sc <- FindNeighbors(sc, dims = 1:40, reduction = "rna.pca",
                    verbose = FALSE)
sc <- FindClusters(sc, resolution = 2,
                   cluster.name = "unintegrated_rna.clusters",
                   verbose = FALSE)
sc <- RunUMAP(sc, dims = 1:40, reduction = "rna.pca",
              reduction.name = "umap.rna.unintegrated",
              verbose = FALSE)
# visualize by batch annotations
p = DimPlot(sc, reduction = "umap.rna.unintegrated", group.by = "patient_id")
ggsave("umap_rna.unintegrated_patient_id.pdf", p, width = 8, height = 6)
ggsave("umap_rna.unintegrated_patient_id.png", p, width = 8, height = 6)

sc <- IntegrateLayers(
  object = sc, method = RPCAIntegration,
  dims = 1:40, assay = "SCT", k.weight = 50,
  normalization.method = "SCT",
  orig.reduction = "rna.pca", new.reduction = "integrated.rna.rpca",
  verbose = FALSE
)
sc = JoinLayers(sc)

sc <- FindNeighbors(sc, reduction = "integrated.rna.rpca",
                    dims = 1:40, verbose = FALSE)
sc <- FindClusters(sc, resolution = 2, cluster.name = "rna.clusters",
                   verbose = FALSE)
sc <- RunUMAP(sc, reduction = "integrated.rna.rpca",
              dims = 1:40, reduction.name = "umap.rna",
              verbose = FALSE)
p <- DimPlot(sc, reduction = "integrated.rna.rpca", group.by = "patient_id")
ggsave("umap_rna.integrated_patient_id.pdf", p, width = 8, height = 6)
ggsave("umap_rna.integrated_patient_id.png", p, width = 8, height = 6)

#### ADT RPCA Batch correction ####
DefaultAssay(sc) = "ADT"
sc[["ADT"]] = split(sc[["ADT"]], f = sc$patient_id)
VariableFeatures(sc) <- rownames(sc[["ADT"]])
sc <- NormalizeData(sc, normalization.method = "CLR", margin = 2,
                    verbose = FALSE)

features <- rownames(sc[["ADT"]])
sc <- ScaleData(sc, features = features,
                do.center = TRUE,
                do.scale = FALSE,
                verbose = FALSE)
sc <- RunPCA(sc, npcs = 40,
             reduction.name = "adt.pca", verbose = FALSE)
sc <- FindNeighbors(sc, dims = 1:40, reduction = "adt.pca",
                    verbose = FALSE)
sc <- FindClusters(sc, resolution = 2,
                   cluster.name = "unintegrated_adt.clusters",
                   verbose = FALSE)
sc <- RunUMAP(sc, dims = 1:40, reduction = "adt.pca",
              reduction.name = "umap.adt.unintegrated",
              verbose = FALSE)
# visualize by batch annotations
p = DimPlot(sc, reduction = "umap.adt.unintegrated", group.by = "patient_id")
ggsave("umap_adt.unintegrated_patient_id.pdf", p, width = 8, height = 6)
ggsave("umap_adt.unintegrated_patient_id.png", p, width = 8, height = 6)

sc <- IntegrateLayers(
  object = sc, method = RPCAIntegration,
  features = features, assay = "ADT", k.weight = 50,
  dims = 1:40, normalization.method = "LogNormalize",
  orig.reduction = "adt.pca", new.reduction = "integrated.adt.rpca",
  verbose = FALSE
)
sc = JoinLayers(sc)

sc <- FindNeighbors(sc, reduction = "integrated.adt.rpca",
                    dims = 1:40, verbose = FALSE)
sc <- FindClusters(sc, resolution = 2, cluster.name = "adt.clusters",
                   verbose = FALSE)
sc <- RunUMAP(sc, reduction = "integrated.adt.rpca",
              dims = 1:40, reduction.name = "umap.adt",
              verbose = FALSE)
p <- DimPlot(sc, reduction = "integrated.adt.rpca", group.by = "patient_id")
ggsave("umap_adt.integrated_patient_id.pdf", p, width = 8, height = 6)
ggsave("umap_adt.integrated_patient_id.png", p, width = 8, height = 6)

#### WNN Integration of RNA and ADT assays ####
if (FALSE) {
  sc <- FindMultiModalNeighbors(
    sc, reduction.list = list("integrated.rna.rpca", "integrated.adt.rpca"),
    dims.list = list(1:40, 1:40), verbose = FALSE
  )
  sc = RunUMAP(sc, nn.name = "weighted.nn",  reduction.name = "wnn.umap",
              reduction.key = "wnnUMAP_", verbose = FALSE)
  p <- DimPlot(sc, reduction = "wnn.umap", group.by = "patient_id")
  ggsave("umap_rna.adt.wnn_patient_id.pdf", p, width = 8, height = 6)
  ggsave("umap_rna.adt.wnn_patient_id.png", p, width = 8, height = 6)
}

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

DefaultAssay(sc) = "RNA"
sc = JoinLayers(sc)
all_markers = FindAllMarkers(sc, verbose = FALSE, assay = "RNA")
all_markers = all_markers[all_markers$p_val_adj < 0.05, ]
write.csv(all_markers, paste("DEG_", graph, ".clusters.res0.25.csv"),
          row.names = FALSE, quote = FALSE)

all_markers = FindAllMarkers(sc, verbose = FALSE, assay = "ADT")
all_markers = all_markers[all_markers$p_val_adj < 0.05, ]
write.csv(all_markers, paste("DEP_", graph, ".clusters.res0.25.csv"),
          row.names = FALSE, quote = FALSE)

non_hc = subset(sc, endotype != "Healthy_Control_Donor")
for (cl_num in unique(Idents(non_hc))) {
  cons_markers = FindConservedMarkers(non_hc, ident.1 = cl_num,
                                      assay = "RNA",
                                      grouping.var = "endotype")
  cons_markers$cluster = cl_num
  if (!exists("all_cons_markers")) {
    all_cons_markers = cons_markers
  } else {
    all_cons_markers = rbind(all_cons_markers, cons_markers)
  }
}
all_cons_markers$gene = rownames(all_cons_markers)
write.csv(all_cons_markers,
          "DEG_nonHC_conserved.markers_cluster.v.endotype.csv",
          quote = FALSE, row.names = FALSE)
###############################################################################
# save Seurat object
saveRDS(sc, paste0(PROJECT_DIR,"/data/qc_rna.hto.adt_", PROJECT_NAME, ".RDS"))
# Save the R session environment information
capture.output(sessionInfo(),
               file=paste0(PROJECT_DIR, "/",
                           PROJECT_NAME,
                           ".Rsession.Info.",
                           gsub("\\D", "", Sys.time()), ".txt"))
###############################################################################
# Optional Additonal Analyses
###############################################################################
# create ShinyCell app with data - MUST pre-authenticate using shinyapps.io
#      token with rsconnect
if (FALSE) {
  adt_cts = LayerData(sc, assay="ADT", layer = "counts")
  sct_cts = LayerData(sc, assay="SCT", layer = "counts")
  adt_data = LayerData(sc, assay="ADT", layer = "data")
  sct_data = LayerData(sc, assay="SCT", layer = "data")
  sct_act_cts = rbind(sct_cts, adt_cts)
  sct_act_data = rbind(sct_data, adt_data)
  sc[["SCT_ADT"]] = CreateAssay5Object(counts = sct_act_cts, data = sct_act_data)
  DefaultAssay(sc) = "SCT_ADT"
  sc = FindVariableFeatures(sc)
  sc_conf = createConfig(sc)
  makeShinyApp(sc, sc_conf, gene.mapping = FALSE,
               shiny.title = paste0(PROJECT_NAME, " RNA + ADT + HTO"),
               shiny.dir = paste0("shiny_",PROJECT_NAME,"_rna"),
               gex.assay="SCT_ADT")
  rsconnect::deployApp(paste0("shiny_",PROJECT_NAME,"_rna"))
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
                        outFile = paste0("qc_sct.adt_",PROJECT_NAME,".h5ad"))
}