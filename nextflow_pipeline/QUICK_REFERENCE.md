# Quick Reference Guide

## Common Commands

### Validate Setup
```bash
# Check all prerequisites
nextflow run main.nf -profile desktop --help

# Run preflight checks only
nextflow run main.nf -profile desktop -entry PREFLIGHT_CHECKS
```

### Desktop Testing
```bash
# Single batch, desktop
./nextflow run main.nf -profile desktop \
    --sample_sheet test_data/sample_sheet.csv \
    --ref_dir /home/Apps/genomes/cellranger \
    --outdir ./results_test

# Resume after failure
./nextflow run main.nf -profile desktop -resume

# Trace execution
./nextflow run main.nf -profile desktop -with-timeline timeline.html
```

### HPC Production
```bash
# Multi-batch SLURM
./nextflow run main.nf -profile slurm \
    --sample_sheet samples_all.csv \
    --ref_dir /mnt/shared/genomes \
    --max_parallel_batches 8 \
    --outdir /scratch/$USER/bcell_atlas

# View live progress
tail -f .nextflow.log

# List running jobs
qstat | grep nextflow  # Or squeue for SLURM
```

### Custom Genomes
```bash
# Mouse (GRCm38)
./nextflow run main.nf --genome GRCm38

# Custom reference (must contain)
# ├── refdata-gex-GRCm38-2024-A/
# ├── refdata-vdj-bcr-GRCm38-alts-ensembl-7.1.0/
# └── feature_refs/
```

## Parameter Cheat Sheet

| Parameter | Default | Notes |
|-----------|---------|-------|
| `sample_sheet` | *required* | CSV with batch_id, FASTQ paths |
| `ref_dir` | `/home/Apps/genomes/cellranger` | GEX + VDJ references |
| `genome` | `GRCh38` | Human (GRCh38/GRCm38/other) |
| `outdir` | `./results` | Output directory |
| `gpu_strategy` | `auto` | auto/manual device selection |
| `include_vdj` | `null` | null=auto-detect, true/false=override |
| `n_hvg` | `2000` | Highly variable genes |
| `harmony_npcs` | `50` | Harmony batch integration (↓ for low mem) |
| `leiden_resolution` | `1.0` | Cluster granularity (↑ = finer clusters) |
| `cellbender_enabled` | `true` | Ambient RNA removal |
| `save_intermediate` | `false` | Save checkpoints (↑ disk, good for debugging) |

## File Structure After Run

```
results/
├── .nextflow/
├── work/                      # Intermediate files (large!)
├── cellranger_outs/
│   ├── batch_1/
│   │   ├── sample_1/
│   │   │   └── outs/
│   │   │       ├── filtered_feature_bc_matrix.h5
│   │   │       ├── web_summary.html
│   │   │       └── metrics_summary.csv
│   │   └── sample_2/outs/
│   └── batch_2/
│       ├── sample_3/outs/
│       └── sample_4/outs/
├── aggregated/
│   ├── batch_1/
│   │   └── aggr_batch_1/outs/
│   │       └── filtered_feature_bc_matrix.h5
│   └── batch_2/aggr_batch_2/outs/
├── qc/
│   ├── cellbender_batch_1.h5ad
│   ├── cellbender_batch_2.h5ad
│   ├── doublet_stats_batch_1.txt
│   └── doublet_stats_batch_2.txt
├── analysis/
│   ├── harmony/
│   │   ├── pre_harmony.h5ad
│   │   └── post_harmony.h5ad
│   ├── clustering/
│   │   ├── umap.h5ad
│   │   └── leiden_clusters.csv
│   └── final/
│       ├── bcell_atlas.h5ad         ← MAIN OUTPUT
│       └── metadata_validation.txt
├── vdj/
│   ├── clonotypes_batch_1.csv
│   ├── clonotypes_batch_2.csv
│   └── vdj_integration_summary.txt
├── metadata/
│   ├── hca_validation.txt
│   └── audit_trail.txt
├── trace.txt                  # Execution timing
├── report.html               # Nextflow report
└── dag.svg                   # Pipeline DAG
```

## Expected Outputs at Each Step

### Cell Ranger (per sample)
```
sample_X/outs/
├── filtered_feature_bc_matrix.h5     ← Used downstream
├── web_summary.html                  ← QC report
├── metrics_summary.csv
├── analysis/
│   └── clustering/
│       └── graphclust/
│           └── clusters.csv
└── logs/
    └── command_log.txt
```

### Aggregated (per batch)
```
aggr_batch_X/outs/
├── filtered_feature_bc_matrix.h5     ← To H5AD
├── aggregation_summary.csv
└── web_summary.html
```

### Final AnnData (.h5ad)
```python
adata.obs.columns
# → ['batch', 'donor_id', 'tissue', 'leiden', 'cellbender_prob', 
#    'is_doublet', 'n_counts', 'n_genes', 'pct_mito', 'clonotype_id', ...]

adata.obsm.keys()
# → ['X_pca', 'X_harmony', 'X_umap']

adata.uns.keys()
# → ['processing', 'hca_compliance', 'vdj_summary', ...]
```

## GPU Troubleshooting

### Check GPU visibility in container
```bash
# Run Singularity with GPU
singularity exec --nv /path/to/rapids-cellranger.sif nvidia-smi

# Output should show GPU(s)
```

### Monitor GPU during run
```bash
# In separate terminal during pipeline execution
watch -n 1 nvidia-smi

# Look for: CUDA processes, memory usage, temperature
```

### If GPU not used
```bash
# 1. Verify GPU availability
nvidia-smi

# 2. Check Singularity GPU binding
export SINGULARITY_BINDPATH="/dev/nvidia*"
export SINGULARITYENV_CUDA_VISIBLE_DEVICES=0,1

# 3. Try manual GPU specification
./nextflow run main.nf --gpu_devices_available 0
```

## Memory Management

### Desktop (limited RAM)
```bash
# Reduce resource usage
./nextflow run main.nf -profile desktop \
    --n_hvg 1000 \          # Fewer genes
    --harmony_npcs 30 \     # Fewer PCs
    --leiden_resolution 0.8 # Coarser clusters
```

### Cleanup intermediate files
```bash
# After successful run, can safely remove:
rm -rf results/work/

# Keep only final outputs and logs
```

## Debugging

### Enable verbose logging
```bash
./nextflow run main.nf -profile desktop \
    -resume \
    -entry CELLRANGER_MULTI_EXECUTE \
    -debug
```

### Re-run single failed process
```bash
# Identify failed process from trace.txt
# E.g., failed task: sha256:abc123...

# Inspect work directory
ls -la results/work/sha256/abc123/

# Re-run with debug
nextflow run main.nf -profile desktop -resume -debug
```

### Inspect intermediate files
```bash
# Load H5AD checkpoint
python3 << 'EOF'
import anndata as ad
adata = ad.read_h5ad("results/analysis/harmony/post_harmony.h5ad")
print(adata)
print("Batches:", adata.obs['batch'].unique())
print("Harmony embedded:", adata.obsm['X_harmony'].shape)
EOF
```

## Performance Optimization

### For large datasets (>1M cells)
```bash
# Enable all GPU features
./nextflow run main.nf \
    --harmony_gpu true \
    --rapids_clustering true \
    --max_parallel_batches 8

# Monitor peak GPU memory
nvidia-smi --query-gpu=index,name,memory.used,memory.total \
    --format=csv --loop=1
```

### For small datasets (<100K cells)
```bash
# Can disable GPU to save container overhead
./nextflow run main.nf \
    --harmony_gpu false \
    --rapids_clustering false
```

### Multi-batch parallelism
```bash
# Desktop (1 batch at a time)
./nextflow run main.nf -profile desktop

# HPC (8 batches in parallel)
./nextflow run main.nf -profile slurm --max_parallel_batches 8
```

## Common Issues & Solutions

| Issue | Cause | Solution |
|-------|-------|----------|
| "reference not found" | Missing genome reference | Check `--ref_dir`, ensure `refdata-gex-*.tar` extracted |
| "GPU out of memory" | Dataset too large for GPU | Reduce `--harmony_npcs`, disable `--rapids_clustering` |
| "cellranger: command not found" | Container issue | Rebuild: `bash scripts/build_containers.sh` |
| "Singularity not found" | Missing Singularity | Install: `apt install singularity-container` |
| Process timeout | Batch too slow | Increase timeout: `conf/desktop.config` maxTime |
| VDJ not included | Auto-detection failed | Force: `--include_vdj true` |

## Next Steps After Completion

### 1. Quality Assessment
```bash
# Load final AnnData
python3 << 'EOF'
import anndata as ad
import scanpy as sc

adata = ad.read_h5ad("results/analysis/final/bcell_atlas.h5ad")

# Basic QC
print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")
print(f"Batches: {adata.obs['batch'].nunique()}")
print(f"Doublets removed: {(adata.obs['is_doublet'] == True).sum()}")
print(f"Clusters: {adata.obs['leiden'].nunique()}")

# Cluster QC
sc.pl.umap(adata, color=['batch', 'leiden', 'donor_id'], ncols=3)
EOF
```

### 2. Downstream Analysis
```r
# Load in Seurat
library(Seurat)
adata <- read_h5ad("results/analysis/final/bcell_atlas.h5ad")
seurat <- as.Seurat(adata)

# Continue with marker detection, annotation, etc.
Idents(seurat) <- "leiden"
markers <- FindAllMarkers(seurat)
```

### 3. Export to Formats
```bash
# Convert to Zarr (large files)
python3 << 'EOF'
import anndata as ad
adata = ad.read_h5ad("bcell_atlas.h5ad")
adata.write_zarr("bcell_atlas.zarr")
EOF

# Export for visualization tools
adata.write_csv("counts.csv")  # For spreadsheets
adata.write_loom("bcell_atlas.loom")  # For Loom viewer
```

## Resources

- **Docs**: See [README.md](README.md)
- **Nextflow**: https://www.nextflow.io/
- **Singularity**: https://sylabs.io/
- **Scanpy**: https://scanpy.readthedocs.io/
- **rapids-singlecell**: https://rapids-singlecell.readthedocs.io/
- **10x Genomics**: https://support.10xgenomics.com/

---

**Save this file for quick reference while running pipelines!**
