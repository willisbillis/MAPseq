# run_phase2.QC_FunF_atac.R
# created by M Elliott Williams (https://github.com/willisbillis) Apr 2024

# Install required packages using the package manager 'pacman'
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)
p_load(Seurat, Signac, GenomeInfoDb, AnnotationHub, biovizBase, ggplot2,
       clustree, dplyr, future, parallel, reticulate, harmony, tidyverse,
       Azimuth, scDblFinder, BiocParallel, aggregation)
p_load_gh("SGDDNB/ShinyCell")
p_load_gh("cellgeni/sceasy")

# Set python path to ensure reticulate packages can be used
python_path = system("which python", intern = TRUE)
use_python(python_path)
# silence random number warning
options(future.rng.onMisuse = "ignore")
set.seed(1234)                # set seed for reproducibility

## Library descriptions ##
# Seurat: functions for single cell data
# Signac: DensityScatter function for QC, ATAC processing
# GenomeInfoDb: database for genomic annotations
# AnnotationHub: database for annotations
# biovizBase: necessary for plots
# ggplot2: functions for plotting
# clustree: plotting clusters vs resolution
# dplyr: pipe command '%>%'
# future: multiprocessing limits
# parallel: multiprocessing limits
# reticulate: set which python to use
# harmony: integration method
# tidyverse: separate function
# Azimuth: cell type annotation
# scDblFinder: doublet detection
# BiocParallel: parallelization for scDblFinder
# aggregation: fisher test for scDblFinder
# ShinyCell: Interact with your data

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
PROJECT_NAME = "FunF_TESTS"
# REPLACE, path to ATAC.ASAP analysis dir from MAPseq pipeline
PROJECT_DIR = paste0("/home/boss_lab/Projects/Scharer_sc/TESTS.MAPseq",
                     "/FunF_TESTS/analysis/ATAC.ASAP")
RAW_SEURAT_PATH = paste0(PROJECT_DIR,
                         "/data/raw_atac.hto_", PROJECT_NAME, ".RDS")
HTO_DEMUX_PATH = paste0(PROJECT_DIR,
                        "/../../pipeline/ATAC.ASAP/hashtag_ref_atac.csv")

GENOME = "GRCh38"                     # REPLACE (GRCh38 or GRCm39)
OUTPUT_FIG_WIDTH =  8                 # inches, width of output figures
OUTPUT_FIG_HEIGHT = 8                 # inches, height of output figures
###############################################################################
#### LOAD DATA ####
###############################################################################
setwd(PROJECT_DIR)
sc_total = readRDS(RAW_SEURAT_PATH)
hto_reference = read.csv(HTO_DEMUX_PATH)
###############################################################################
#### ATTACH LATEST GENE ANNOTATIONS TO ATAC DATA ####
###############################################################################
DefaultAssay(sc_total) = "ATAC"
# Add the gene annotations to the peaks in the assay
ah <- AnnotationHub()
if (GENOME %in% c("GRCm39")) {
  ahDb <- query(ah, pattern = c("Mus Musculus", "EnsDb"))
} else if (GENOME %in% c("GRCh38")) {
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
annotations = GetGRangesFromEnsDb(ensdb = GENOME_DATA)

# change to UCSC style
annotations <- renameSeqlevels(annotations,
                               mapSeqlevels(seqlevels(annotations), "UCSC"))

# add the gene information to the object
Annotation(sc_total) = annotations
###############################################################################
#### CALCULATE QC METRICS (ATAC) ####
###############################################################################
# compute nucleosome signal score per cell
sc_total = NucleosomeSignal(object = sc_total, verbose = FALSE)

# compute TSS enrichment score per cell
sc_total = TSSEnrichment(object = sc_total, verbose = FALSE)

# add blacklist ratio and fraction of reads in peaks
sc_total$pct_reads_in_peaks = sc_total$peak_region_fragments /
  sc_total$passed_filters * 100
sc_total$blacklist_ratio = sc_total$blacklist_region_fragments /
  sc_total$peak_region_fragments

# Doublet Detection
sce = scDblFinder(as.SingleCellExperiment(sc_total), artificialDoublets = 1,
                  aggregateFeatures = TRUE, samples = "atac_id",
                  nfeatures = 25, processing = "normFeatures",
                  BPPARAM = MulticoreParam(max_cores))

to_exclude <- GRanges(c("M", "chrM", "MT", "X", "Y", "chrX", "chrY"),
                      IRanges(1L, width = 10^8))
res <- amulet(Fragments(sc_total)[[1]]@path, regionsToExclude = to_exclude)
res$scDblFinder.p <- 1 - colData(sce)[row.names(res), "scDblFinder.score"]
res$combined <- apply(res[, c("scDblFinder.p", "p.value")], 1,
  FUN = function(x) {
    x[x < 0.001] <- 0.001 # prevent too much skew from very small or 0 p-values
    suppressWarnings(aggregation::fisher(x))
  }
)
sc_total$scDblFinder.score <- res$combined

p = DensityScatter(sc_total, "peak_region_fragments", "scDblFinder.score",
                   quantiles = TRUE, log_x = TRUE, log_y = FALSE)
ggsave("scatter_peakfrags.v.dbl_alldata.png",
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)

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
#### ATAC QC CUTOFFS ####
###############################################################################
# PAUSE, view scatter figures above and determine appropriate cutoffs below
MIN_PEAK_FRAGMENTS = 500   # REPLACE, minimum peak fragments per cell
MIN_PCT_RiP = 65            # REPLACE, minimum percent reads in peaks per cell
MAX_BLACKLIST_RATIO = 1.0   # REPLACE, maximum blacklist ratio per cell
MAX_NUCLEOSOME_SIG = 1      # REPLACE, maximum nucleosome signal per cell
MIN_TSS = 4                 # REPLACE, minimum TSS enrichment score per cell
DBL_LIMIT = 0.25            # REPLACE, minimum Doublet score accepted

p = DensityScatter(sc_total, "peak_region_fragments", "TSS.enrichment",
                   quantiles = TRUE, log_x = TRUE, log_y = FALSE)
p = p +
  geom_hline(yintercept = MIN_TSS, linetype = "dashed") +
  geom_vline(xintercept = MIN_PEAK_FRAGMENTS,
             linetype = "dashed")
ggsave("scatter_peakfrags.v.TSSe_filtered.png",
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)
###############################################################################
#### QUANTIFY QC FILTERING ####
###############################################################################
# adjust metadata to accomodate Seurat's AggregateExpression
sc_total$atac_id = gsub("_", "-", sc_total$atac_id)
sc_total$patient_id = gsub("_", "-", sc_total$patient_id)
hto_reference$library_id = gsub("_", "-", hto_reference$library_id)
hto_reference$patient_id = gsub("_", "-", hto_reference$patient_id)
# create new column for unique sample ID - adjust as needed for each dataset
hto_reference$sample_id = paste(hto_reference$library_id,
                                hto_reference$patient_id,
                                sep = "-")
sc_total$sample_id = paste(sc_total$atac_id,
                           sc_total$patient_id,
                           sep = "-")
stats = data.frame(sample_id = hto_reference$sample_id)
sample_id_counts = as.data.frame(table(sc_total$sample_id))
stats$Unfiltered_Cells = sample_id_counts$Freq[match(stats$sample_id,
                                                     sample_id_counts$Var1)]
stats[is.na(stats)] = 0
cell_ct = AggregateExpression(sc_total,
                              group.by = "sample_id")$ATAC %>%
  colSums
stats$Unfiltered_Avg_Accessibility = round(cell_ct[match(stats$sample_id,
                                                         names(cell_ct))] /
                                             stats$Unfiltered_Cells)
cell_ct = AggregateExpression(sc_total,
                              group.by = "sample_id")$HTO %>%
  colSums
stats$Unfiltered_Avg_HTO = round(cell_ct[match(stats$sample_id,
                                               names(cell_ct))] /
                                   stats$Unfiltered_Cells)

sc <- subset(sc_total,
             peak_region_fragments > MIN_PEAK_FRAGMENTS &
               pct_reads_in_peaks > MIN_PCT_RiP &
               blacklist_ratio < MAX_BLACKLIST_RATIO &
               nucleosome_signal < MAX_NUCLEOSOME_SIG &
               TSS.enrichment > MIN_TSS &
               scDblFinder.score > DBL_LIMIT)

sample_id_counts = as.data.frame(table(sc$sample_id))
stats$Filtered_Cells = sample_id_counts$Freq[match(stats$sample_id,
                                                   sample_id_counts$Var1)]
stats[is.na(stats)] = 0

cell_ct = AggregateExpression(sc,
                              group.by = "sample_id")$ATAC %>%
  colSums
stats$Filtered_Avg_Accessibility = round(cell_ct[match(stats$sample_id,
                                                       names(cell_ct))] /
                                           stats$Filtered_Cells)
cell_ct = AggregateExpression(sc,
                              group.by = "sample_id")$HTO %>%
  colSums
stats$Filtered_Avg_HTO = round(cell_ct[match(stats$sample_id,
                                             names(cell_ct))] /
                                 stats$Filtered_Cells)
stats[is.na(stats)] = 0

write.csv(stats, "QC.atac_sampleID_filtering.stats.csv",
          quote = FALSE, row.names = FALSE)
###############################################################################
#### SAVE RAW SEURAT OBJECT ####
###############################################################################
saveRDS(sc_total,
        paste0(PROJECT_DIR, "/data/raw_atac.hto_", PROJECT_NAME, ".RDS"))
###############################################################################
#### CREATE DEMONSTRATION FIGURES ####
p = FeatureScatter(sc_total, "nCount_HTO", "nCount_ATAC", group.by = "asap_id",
                   split.by = "patient_id") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  stat_ellipse(aes(group = asap_id), sc_total@meta.data)
ggsave("featscatter_nctHTO.v.nctATAC_patient_id_alldata.png", p,
       width = 10, height = 4)
p = FeatureScatter(sc, "nCount_HTO", "nCount_ATAC", group.by = "asap_id",
                   split.by = "patient_id") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  stat_ellipse(aes(group = asap_id), sc_total@meta.data)
ggsave("featscatter_nctHTO.v.nctATAC_patient_id_filtereddata.png", p,
       width = 10, height = 4)

cells_use = colnames(sc_total[, !is.na(sc_total$HTO_maxID)])
p = FeatureScatter(sc_total, "nCount_HTO", "nCount_ATAC", group.by = "asap_id",
                   split.by = "HTO_maxID", cells = cells_use, ncol = ncol) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  stat_ellipse(aes(group = asap_id),
               sc_total@meta.data[!is.na(sc_total@meta.data$HTO_maxID), ])
ggsave("featscatter_nctHTO.v.nctATAC_HTO_maxid_alldata.png", p,
       width = OUTPUT_FIG_WIDTH * 2, height = OUTPUT_FIG_HEIGHT)

cells_use = colnames(sc[, !is.na(sc$HTO_maxID)])
p = FeatureScatter(sc, "nCount_HTO", "nCount_ATAC", group.by = "asap_id",
                   split.by = "HTO_maxID", cells = cells_use, ncol = ncol) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  stat_ellipse(aes(group = asap_id),
               sc@meta.data[!is.na(sc@meta.data$HTO_maxID), ])
ggsave("featscatter_nctHTO.v.nctATAC_HTO_maxid_filtereddata.png", p,
       width = OUTPUT_FIG_WIDTH * 2, height = OUTPUT_FIG_HEIGHT)
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