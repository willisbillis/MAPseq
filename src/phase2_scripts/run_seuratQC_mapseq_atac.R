# run_seuratQC_mapseq_atac.R - created by M Elliott Williams (https://github.com/willisbillis) Feb 2024

# Install required packages using the package manager 'pacman'
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)
p_load(Seurat, Signac, GenomeInfoDb, AnnotationHub, biovizBase, ggplot2,
       clustree, dplyr, future, parallel, reticulate, harmony, multtest,
       metap)
p_load_gh("SGDDNB/ShinyCell")
p_load_gh("cellgeni/sceasy")

# Set python path to ensure reticulate packages can be used
python_path = system("which python", intern = TRUE)
use_python(python_path)
# make sure you are running Seurat v5
options(Seurat.object.assay.version = "v5")
# silence random number warning
options(future.rng.onMisuse = "ignore")
set.seed(1234)                # set seed for reproducibility
## Library descriptions ##
# Seurat: functions for single cell data
# Signac: DensityScatter function for QC
# GenomeInfoDb: database for genomic annotations
# EnsDb.Mmusculus.v79: database for mm10 annotations
# EnsDb.Hsapiens.v86: database for hg38 annotations
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
# REPLACE, path to ATAC.ASAP analysis dir from MAPseq pipeline
PROJECT_DIR = "/home/boss_lab/Projects/Scharer_sc/LRA.MAPseq/LRA_all/analysis/ATAC.ASAP"
RAW_SEURAT_PATH = paste0(PROJECT_DIR,"/data/raw_atac.hto_", PROJECT_NAME, ".RDS")
HTO_DEMUX_PATH = paste0(PROJECT_DIR,"/../../pipeline/ATAC.ASAP/hashtag_ref_atac.csv")

GENOME = "GRCh38"                     # REPLACE (GRCh38, GRCm38, or GRCm39)
OUTPUT_FIG_WIDTH =  8               # inches, width of output figures
OUTPUT_FIG_HEIGHT = 8

###############################################################################
setwd(PROJECT_DIR)
sc_total = readRDS(RAW_SEURAT_PATH)
###############################################################################
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

top_confusion_matrix = as.data.frame(table(sc_total$HTO_maxID[sc_total$HTO_margin > quantile(sc_total$HTO_margin, 0.95, na.rm=T)],
       sc_total$HTO_secondID[sc_total$HTO_margin > quantile(sc_total$HTO_margin, 0.95, na.rm=T)]))
bot_confusion_matrix = as.data.frame(table(sc_total$HTO_maxID[sc_total$HTO_margin < quantile(sc_total$HTO_margin, 0.05, na.rm=T)],
       sc_total$HTO_secondID[sc_total$HTO_margin < quantile(sc_total$HTO_margin, 0.05, na.rm=T)]))

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
write.csv(confusion_matrix_all, paste0("HTB.combos_",PROJECT_NAME,"_metrics.csv"), quote = F, row.names = F)

Idents(sc_total) = "hash.ID"

p = FeatureScatter(sc_total, cells = colnames(sc_total)[sc_total$hash.ID %in% c(best_htos, "Doublet")],
    feature1 = best_htos[1], feature2 = best_htos[2])
ggsave(paste0("scatter_best.hto.separation_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT, width = OUTPUT_FIG_WIDTH)

p = FeatureScatter(sc_total, cells = colnames(sc_total)[sc_total$hash.ID %in% c(worst_htos, "Doublet")],
    feature1 = worst_htos[1], feature2 = worst_htos[2])
ggsave(paste0("scatter_worst.hto.separation_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT, width = OUTPUT_FIG_WIDTH)
###############################################################################
DefaultAssay(sc_total) = "ATAC"
# Add the gene annotations to the peaks in the assay
ah <- AnnotationHub()
if (GENOME %in% c("GRCm38", "GRCm39")) {
  ahDb <- query(ah, pattern = c("Mus Musculus", "EnsDb"))
} else if (GENOME %in% c("GRCh38", "GRCh39")) {
  ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb"))
}  else {
  print(paste("[ERROR] Requested genome,",
              GENOME,
              ",not available or not implemented. Exiting."))
  quit()
}

ahDb = ahDb[ahDb$genome == GENOME] # Filter to specific genome
ahDb = ahDb[ahDb$rdatadateadded == ahDb@date] # Find latest version
GENOME_DATA = ahDb[[names(ahDb@.db_uid)]]

# extract gene annotations from EnsDb
annotations = GetGRangesFromEnsDb(ensdb=GENOME_DATA)

# change to UCSC style
annotations <- renameSeqlevels(annotations, mapSeqlevels(seqlevels(annotations), "UCSC"))

# add the gene information to the object
Annotation(sc_total) = annotations
###############################################################################
# Compute some basic QC metrics on the assay

# compute nucleosome signal score per cell
sc_total = NucleosomeSignal(object=sc_total, verbose = FALSE)

# compute TSS enrichment score per cell
sc_total = TSSEnrichment(object=sc_total, verbose = FALSE)

# add blacklist ratio and fraction of reads in peaks
sc_total$pct_reads_in_peaks = sc_total$peak_region_fragments /
  sc_total$passed_filters * 100
sc_total$blacklist_ratio = sc_total$blacklist_region_fragments /
  sc_total$peak_region_fragments

p = DensityScatter(sc_total, "peak_region_fragments", "TSS.enrichment",
                   quantiles = TRUE, log_x = TRUE, log_y = FALSE)
ggsave("scatter_peakfrags.v.TSSe_alldata.png",
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)

p = DensityScatter(sc_total, "peak_region_fragments", "pct_reads_in_peaks",
                   quantiles = TRUE, log_x = TRUE, log_y = FALSE)
ggsave("scatter_peakfrags.v.prip_alldata.png",
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)

p = DensityScatter(sc_total, "peak_region_fragments", "nucleosome_signal",
                   quantiles = TRUE, log_x = TRUE, log_y = FALSE)
ggsave("scatter_peakfrags.v.nsig_alldata.png",
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)

p = DensityScatter(sc_total, "peak_region_fragments", "blacklist_ratio",
                   quantiles = TRUE, log_x = TRUE, log_y = FALSE)
ggsave("scatter_peakfrags.v.blacklist_alldata.png",
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)
###############################################################################
###############################################################################
# PAUSE, view scatter figures above and determine appropriate cutoffs below
MIN_PEAK_FRAGMENTS = 1000   # REPLACE, minimum peak fragments per cell
MAX_PEAK_FRAGMENTS = 35000  # REPLACE, maximum peak fragments per cell
MIN_PCT_RiP = 80            # REPLACE, minimum percent reads in peaks
MAX_BLACKLIST_RATIO = 1.0   # REPLACE, maximum blacklist ratio per cell
MAX_NUCLEOSOME_SIG = 1      # REPLACE, maximum nucleosome signal per cell
MIN_TSS = 2                 # REPLACE, minimum TSS enrichment score per cell

p = DensityScatter(sc_total, "peak_region_fragments", "TSS.enrichment",
                   quantiles = TRUE, log_x = TRUE, log_y = FALSE)
p = p +
  geom_hline(yintercept=MIN_TSS, linetype = "dashed") +
  geom_vline(xintercept=c(MIN_PEAK_FRAGMENTS, MAX_PEAK_FRAGMENTS),
             linetype = "dashed")
ggsave("scatter_peakfrags.v.TSSe_filtered.png",
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)
###############################################################################
###############################################################################
# Filtering of cells based on QC criteria

# Expression Filter
hto_reference = read.csv(HTO_DEMUX_PATH)
stats = data.frame(Patients = hto_reference$patient_id)
patient_id.counts = as.data.frame(table(sc_total$patient_id))
stats$Unfiltered_Cells = patient_id.counts$Freq[match(stats$Patients, patient_id.counts$Var1)]
stats[is.na(stats)] = 0
cell_counts = round(AggregateExpression(sc_total, group.by="patient_id")$ATAC %>% colSums)
stats$Unfiltered_Average_Accessibility = cell_counts[match(stats$Patients, names(cell_counts))] /
       stats$Unfiltered_Cells

sc <- subset(sc_total,
             peak_region_fragments > MIN_PEAK_FRAGMENTS &
               peak_region_fragments < MAX_PEAK_FRAGMENTS &
               pct_reads_in_peaks > MIN_PCT_RiP &
               blacklist_ratio < MAX_BLACKLIST_RATIO &
               nucleosome_signal < MAX_NUCLEOSOME_SIG &
               TSS.enrichment > MIN_TSS)

# HTO Filter
sc <- subset(sc, subset = HTO_classification.global != "Doublet")

patient_id.counts = as.data.frame(table(sc$patient_id))
stats$Filtered_Cells = patient_id.counts$Freq[match(stats$Patients, patient_id.counts$Var1)]
stats[is.na(stats)] = 0
cell_counts = round(AggregateExpression(sc, group.by="patient_id")$ATAC %>% colSums)
stats$Filtered_Average_Accessibility = cell_counts[match(stats$Patients, names(cell_counts))] /
       stats$Filtered_Cells
stats[is.na(stats)] = 0

write.csv(stats, "QC.atac_patientID_filtering.stats.csv", quote = F, row.names = F)
###############################################################################
# Clustering, DAR (differentially expressed genes) for ATAC assay
DefaultAssay(sc) = "ATAC"
batch_column = "endotype"
n_dims = 30

sc_na = sc[, colnames(sc)[is.na(sc$endotype)]]
sc <- sc[, colnames(sc)[!is.na(sc$endotype)]]
sc_na$label_exists = FALSE
sc$label_exists = TRUE

sc <- RunTFIDF(sc, min.cells = 1)
sc <- FindTopFeatures(sc, min.cutoff = "q0")
sc <- RunSVD(sc, n = n_dims, reduction.name = "atac.lsi",
             verbose = FALSE)

sc <- FindNeighbors(sc, dims = 2:n_dims, reduction = "atac.lsi",
                    verbose = FALSE)
sc <- RunUMAP(sc, dims = 2:n_dims, reduction = "atac.lsi",
              reduction.name = "umap.atac.unintegrated",
              return.model = TRUE, verbose = FALSE)
# visualize by batch annotations
p = DimPlot(sc, reduction = "umap.atac.unintegrated", group.by = batch_column)
ggsave(paste0("umap_atac.unintegrated_", batch_column, ".pdf"), p,
       width = 8, height = 6)
ggsave(paste0("umap_atac.unintegrated_", batch_column, ".png"), p,
       width = 8, height = 6)

sc <- RunHarmony(
  object = sc,
  group.by.vars = batch_column,
  reduction.use = "atac.lsi",
  assay.use = "ATAC",
  project.dim = FALSE
)
sc <- FindNeighbors(sc, reduction = "harmony",
                    dims = 2:n_dims, verbose = FALSE)
sc <- RunUMAP(sc, reduction = "harmony",
              dims = 2:n_dims, reduction.name = "umap.atac",
              return.model = TRUE, verbose = FALSE)

p <- DimPlot(sc, reduction = "umap.atac", group.by = batch_column)
ggsave(paste0("umap_atac.integrated_", batch_column, "_known.labels.pdf"), p,
       width = 8, height = 6)
ggsave(paste0("umap_atac.integrated_", batch_column, "_known.labels.png"), p,
       width = 8, height = 6)

# find transfer anchors
sc_na = RunTFIDF(sc_na, min.cells = 1)
sc_na = FindTopFeatures(sc_na, min.cutoff = "q0")
sc_na <- RunSVD(sc_na, n = n_dims, reduction.name = "atac.lsi",
                verbose = FALSE)

transfer_anchors <- FindTransferAnchors(
  reference = sc,
  query = sc_na,
  reference.reduction = "atac.lsi",
  reduction = "lsiproject",
  dims = 2:n_dims
)

# map query onto the reference dataset
sc_na <- MapQuery(
  anchorset = transfer_anchors,
  reference = sc,
  query = sc_na,
  refdata = sc$endotype,
  reference.reduction = "atac.lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = "umap.atac.unintegrated"
)
sc_na[[batch_column]] = sc_na$predicted.id
sc <- merge(sc, sc_na)

# Run dimensionality reduction and clustering on the new
#      fully annotated dataset
sc <- RunTFIDF(sc, min.cells = 1)
sc <- FindTopFeatures(sc, min.cutoff = "q0")
sc <- RunSVD(sc, n = n_dims, reduction.name = "atac.lsi",
             verbose = FALSE)
sc <- FindNeighbors(sc, dims = 2:n_dims, reduction = "atac.lsi",
                    verbose = FALSE)
sc <- FindClusters(sc, resolution = 2,
                   cluster.name = "unintegrated_atac.clusters",
                   verbose = FALSE)
sc <- RunUMAP(sc, dims = 2:n_dims, reduction = "atac.lsi",
              reduction.name = "umap.atac.unintegrated",
              verbose = FALSE)

# visualize by batch annotations
p = DimPlot(sc, reduction = "umap.atac.unintegrated", group.by = batch_column,
            split.by = "label_exists")
ggsave(paste0("umap_atac.unintegrated_", batch_column, "_label.split.pdf"), p,
       width = 16, height = 6)
ggsave(paste0("umap_atac.unintegrated_", batch_column, "_label.split.png"), p,
       width = 16, height = 6)

sc <- RunHarmony(
  object = sc,
  group.by.vars = batch_column,
  reduction.use = "atac.lsi",
  assay.use = "ATAC",
  project.dim = FALSE
)
sc <- FindNeighbors(sc, reduction = "harmony",
                    dims = 2:n_dims, verbose = FALSE)
sc <- RunUMAP(sc, reduction = "harmony",
              dims = 2:n_dims, reduction.name = "umap.atac",
              return.model = TRUE, verbose = FALSE)

p <- DimPlot(sc, reduction = "umap.atac", group.by = batch_column,
             split.by = "label_exists")
ggsave(paste0("umap_atac.integrated_", batch_column, ".pdf"), p,
       width = 16, height = 6)
ggsave(paste0("umap_atac.integrated_", batch_column, ".png"), p,
       width = 16, height = 6)

# Different cluster resolutions for ATAC
graph = "ATAC_snn"
for (res in c(1, 0.5, 0.25, 0.1, 0.05)) {
  sc <- FindClusters(sc, resolution = res, graph.name = graph,
                    algorithm = 3, verbose = FALSE)
}

p = clustree(sc, prefix = paste0(graph, "_res."))
ggsave("clustree_clusters_rna.png", p,
       width=OUTPUT_FIG_WIDTH, height=OUTPUT_FIG_HEIGHT)

sc$seurat_clusters = sc[[paste0(graph, "_res.", 0.25)]]
Idents(sc) = "seurat_clusters"
sc$seurat_clusters = factor(sc$seurat_clusters)

all_markers = FindAllMarkers(sc, verbose = FALSE, assay = "ATAC")
all_markers = all_markers[all_markers$p_val_adj < 0.05, ]
closest_feats = ClosestFeature(sc, regions=rownames(all_markers))
all_markers$gene = closest_feats$gene_name[match(rownames(all_markers),
                                                 closest_feats$query_region)]
write.csv(all_markers, paste("DAR_", graph, ".clusters.res0.25.csv"),
          row.names = FALSE, quote = FALSE)

non_hc = subset(sc, endotype != "Healthy_Control_Donor")
for (cl_num in unique(Idents(sc))) {
  cons_markers = FindConservedMarkers(non_hc, ident.1 = cl_num,
                                      assay = "ATAC",
                                      grouping.var = "endotype")
  closest_feats = ClosestFeature(non_hc, regions=rownames(cons_markers))
  cons_markers$gene = closest_feats$gene_name[match(rownames(cons_markers),
                                                    closest_feats$query_region)]
  cons_markers$cluster = cl_num
  if (!exists("all_cons_markers")) {
    all_cons_markers = cons_markers
  } else {
    all_cons_markers = rbind(all_cons_markers, cons_markers)
  }
}
write.csv(all_cons_markers,
          "DAR_nonHC_conserved.markers_cluster.v.endotype.csv",
          quote = FALSE, row.names = FALSE)
###############################################################################
# save Seurat object
saveRDS(sc, paste0(PROJECT_DIR,"/data/qc_atac.hto_", PROJECT_NAME, ".RDS"))
# Save the R session environment information
capture.output(sessionInfo(),
               file=paste0(PROJECT_DIR, "/",
                           PROJECT_NAME,
                           ".Rsession.Info.",
                           gsub("\\D", "", Sys.time()), ".txt"))
###############################################################################
###############################################################################
# create ShinyCell app with data - MUST pre-authenticate using shinyapps.io
#      token with rsconnect
if (FALSE) {
  sc <- AddMetaData(sc, t(LayerData(sc, assay="HTO")))
  DefaultAssay(sc) = "ATAC"
  gene_activities = GeneActivity(sc)
  sc[["pseudoRNA"]] <- CreateAssayObject(counts = gene_activities)
  sc <- NormalizeData(
    object = sc,
    assay = "pseudoRNA",
    normalization.method = "LogNormalize",
    scale.factor = median(sc$nCount_pseudoRNA)
  )
  DefaultAssay(sc) = "pseudoRNA"
  sc <- FindVariableFeatures(sc)
  sc_conf = createConfig(sc)
  makeShinyApp(sc, sc_conf, gene.mapping = TRUE,
               shiny.title = paste0(PROJECT_NAME, " ATAC (pseudoRNA) + HTO"),
               shiny.dir = paste0("shiny_", PROJECT_NAME, "_atac"),
               gex.assay="pseudoRNA")
  rsconnect::deployApp(paste0("shiny_", PROJECT_NAME, "_atac"))
}