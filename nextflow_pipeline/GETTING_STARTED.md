# Deployment Checklist & Getting Started

## ✅ Pre-Deployment Checklist

### System Requirements
- [ ] Linux system (Ubuntu 20.04+, CentOS 7+, or RedHat equivalent)
- [ ] Java 11+ installed (`java -version`)
- [ ] Python 3.8+ installed (`python3 --version`)
- [ ] GPU support (optional but recommended):
  - [ ] NVIDIA GPU with Compute Capability 7.0+ (V100, A100, RTX 30xx/40xx series)
  - [ ] NVIDIA drivers installed (`nvidia-smi` works)
  - [ ] CUDA 11.0+ runtime
  - [ ] cuDNN installed

### Software Installation
- [ ] Nextflow installed (`nextflow -version`)
- [ ] Singularity 3.7+ installed (`singularity --version`)
- [ ] Git installed (`git --version`)
- [ ] Cell Ranger 8.0+ installed (`cellranger --version`)
- [ ] Python packages:
  ```bash
  pip install anndata scanpy pandas numpy scipy matplotlib seaborn
  ```

### File System Preparation
- [ ] Reference directory exists with Cell Ranger references:
  - [ ] `refdata-gex-GRCh38-2024-A/` (or custom genome)
  - [ ] `refdata-vdj-bcr-GRCh38-alts-ensembl-7.1.0/` (for VDJ)
  - [ ] Feature reference CSV (e.g., HTOB_feature_ref.csv)
- [ ] FASTQ input directory prepared with all sequencing data
- [ ] Output directory exists and has write permissions
- [ ] Singularity cache directory accessible:
  - Desktop: `~/.singularity/cache` (must have 10+ GB free)
  - HPC: `/mnt/shared/singularity/cache` (shared storage)

## 📋 First-Time Setup (5-10 minutes)

### Step 1: Clone Repository
```bash
cd /home/elliott/Apps/github/MAPseq
git clone https://github.com/yourusername/bcell-atlas-pipeline nextflow_pipeline
cd nextflow_pipeline
chmod +x nextflow
```

### Step 2: Verify Installation
```bash
# Check Nextflow
./nextflow -version

# Check Singularity
singularity --version

# Check GPU (if available)
nvidia-smi

# Verify directory structure
ls -la
# Expected: main.nf, nextflow.config, params.config, modules/, conf/, containers/, scripts/, test_data/
```

### Step 3: Build Containers (30-60 minutes)
```bash
# Build GPU-accelerated containers
bash scripts/build_containers.sh ~/.singularity/cache containers/

# Verify builds
ls -lh ~/.singularity/cache/*.sif

# Test a container
singularity exec --nv ~/.singularity/cache/rapids-cellranger.sif nvidia-smi
```

### Step 4: Prepare Sample Sheet
```bash
# Create samples.csv in project root with your data:
# sample_id,batch_id,gex_fastq_r1,gex_fastq_r2,...

# Or copy test data
cp test_data/sample_sheet.csv my_samples.csv
# Edit with your actual sample data
```

### Step 5: Configure References (Desktop Example)
```bash
# Check your reference location
ls /home/Apps/genomes/cellranger/

# Expected structure:
# ├── refdata-gex-GRCh38-2024-A/
# ├── refdata-vdj-bcr-GRCh38-alts-ensembl-7.1.0/
# └── features/
#     ├── HTOB_feature_ref.csv
#     └── TSC_feature_ref.csv
```

## 🚀 Running Your First Pipeline (Desktop)

### Quick Start (Test Run)
```bash
# Run test with provided sample data (~30 min on desktop)
./nextflow run main.nf -profile desktop \
    --sample_sheet test_data/sample_sheet.csv \
    --ref_dir /home/Apps/genomes/cellranger \
    --outdir ./results_test \
    -with-timeline timeline.html

# Monitor progress
tail -f .nextflow.log

# View results
ls -la results_test/
```

### Production Run (Your Data)
```bash
# Prepare your samples.csv with actual paths

./nextflow run main.nf -profile desktop \
    --sample_sheet my_samples.csv \
    --ref_dir /home/Apps/genomes/cellranger \
    --outdir ./results \
    --n_hvg 3000 \
    --harmony_npcs 50 \
    -with-timeline timeline.html \
    -with-report report.html
```

### Monitor Execution
```bash
# In separate terminal, monitor GPU (if using GPU)
watch -n 1 nvidia-smi

# Monitor memory
watch -n 1 free -h

# View pipeline DAG
# Results: results/dag.svg (open in browser after run completes)
```

## 💻 HPC Setup (SLURM)

### Cluster-Specific Configuration
```bash
# Edit conf/hpc.config for your cluster:
process {
    executor = 'slurm'
    queue = 'gpu-compute'      # Your GPU queue name
    time = { 8.h * task.attempt }
    memory = { 64.GB * task.attempt }
    cpus = 16
    
    withLabel: 'gpu_*' {
        clusterOptions = '--gres=gpu:1'  # Adjust for your cluster
    }
}
```

### Submit HPC Job
```bash
# Create job submission script: submit.sh
#!/bin/bash
#SBATCH --job-name=bcell_atlas
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --gpus=4
#SBATCH --mem=256G
#SBATCH --output=logs/bcell_%J.log

cd /scratch/$USER/bcell_atlas_work

./nextflow run /path/to/pipeline/main.nf \
    -profile slurm \
    --sample_sheet /path/to/samples.csv \
    --ref_dir /mnt/shared/genomes \
    --outdir /scratch/$USER/results \
    --max_parallel_batches 8 \
    -resume

# Submit
sbatch submit.sh

# Monitor
squeue -u $USER
tail -f logs/bcell_*.log
```

## 📊 Expected Output Structure

After successful completion:
```
results/
├── cellranger_outs/              # Per-sample Cell Ranger outputs
├── aggregated/                   # Batch-aggregated matrices
├── qc/                          # Quality control reports
│   ├── cellbender/
│   └── doublets/
├── analysis/
│   ├── harmony/
│   ├── clustering/
│   └── final/
│       └── bcell_atlas.h5ad     ← MAIN OUTPUT
├── vdj/                         # VDJ clonotype analysis (if included)
├── metadata/                    # HCA compliance reports
├── trace.txt                    # Execution metrics
├── report.html                  # Nextflow report
└── dag.svg                      # Pipeline DAG
```

**Load main output:**
```python
import anndata as ad
adata = ad.read_h5ad("results/analysis/final/bcell_atlas.h5ad")
print(adata)  # Should show 1M+ cells × 36K genes
```

## 🔧 Common Setup Issues

### Issue: Container build fails
```bash
# Solution 1: Check internet connectivity
ping nvidia.com

# Solution 2: Clear cache and retry
rm -rf ~/.singularity/cache/*.sif
bash scripts/build_containers.sh ~/.singularity/cache containers/

# Solution 3: Use pre-built containers (if available)
# Contact team for shared container location
```

### Issue: GPU not detected in Nextflow
```bash
# Solution: Set environment variables
export SINGULARITY_BINDPATH="/dev/nvidia*"
export SINGULARITYENV_CUDA_VISIBLE_DEVICES=0,1
./nextflow run main.nf ...
```

### Issue: Reference directory not found
```bash
# Solution: Verify path
ls -la /home/Apps/genomes/cellranger/
# If missing, download or create references:
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
tar -xzf refdata-gex-GRCh38-2024-A.tar.gz -C /home/Apps/genomes/cellranger/
```

### Issue: Out of disk space
```bash
# Check available space
df -h /scratch  # or your output directory

# Solution: Clean intermediate files after successful runs
rm -rf results/work/
# (Keep trace.txt, report.html, final results)
```

### Issue: Memory errors
```bash
# For desktop with limited RAM:
./nextflow run main.nf -profile desktop \
    --n_hvg 1500 \           # Reduce features
    --harmony_npcs 30 \      # Reduce PCs
    --save_intermediate false
```

## 📈 Next: Processing Your Data

### Sample Sheet Format
Create `my_samples.csv` with columns:
```csv
sample_id,batch_id,gex_fastq_r1,gex_fastq_r2,gex_feature_bc,hto_fastq_r1,hto_fastq_r2,vdj_fastq_r1,vdj_fastq_r2,donor_id,tissue,facs_gates
sample_1,batch_1,/data/sample_1_R1.fq.gz,/data/sample_1_R2.fq.gz,HTOB_feature_ref.csv,/data/sample_1_HTO_R1.fq.gz,/data/sample_1_HTO_R2.fq.gz,/data/sample_1_VDJ_R1.fq.gz,/data/sample_1_VDJ_R2.fq.gz,donor_A,blood,CD19+CD27-
```

**Column guide:**
- `sample_id`: Unique identifier
- `batch_id`: Batch grouping (samples in same batch will be aggregated together)
- `gex_fastq_r1/r2`: Gene expression FASTQ files (R1=read1, R2=read2/barcode)
- `gex_feature_bc`: Feature barcode reference CSV (for HTO detection)
- `hto_fastq_r1/r2`: HTO (sample hashing) FASTQ (optional, leave empty if not used)
- `vdj_fastq_r1/r2`: VDJ (TCR/BCR) FASTQ (optional, auto-skip if empty)
- `donor_id`: Donor identifier (for merging or tracking)
- `tissue`: Tissue source (e.g., blood, spleen, lymph_node)
- `facs_gates`: FACS sorting criteria applied (e.g., "CD19+CD27-")

### Run Production Pipeline
```bash
./nextflow run main.nf -profile desktop \
    --sample_sheet my_samples.csv \
    --ref_dir /home/Apps/genomes/cellranger \
    --outdir ./results \
    --n_hvg 3000 \
    --harmony_npcs 50 \
    --leiden_resolution 1.0 \
    --cellbender_enabled true \
    --include_vdj auto \
    --save_intermediate true
```

## ✨ Advanced Customization

### Custom Genome
```bash
# Mouse
./nextflow run main.nf --genome GRCm38

# Custom references
./nextflow run main.nf --ref_dir /custom/refs
```

### GPU Optimization
```bash
# Manual GPU specification
./nextflow run main.nf --gpu_devices_available 0,1

# Disable GPU (for testing CPU)
./nextflow run main.nf --harmony_gpu false
```

### Output Customization
```bash
# Save all intermediate checkpoints
./nextflow run main.nf --save_intermediate true

# This creates: pre_harmony.h5ad, post_harmony.h5ad, etc.
```

## 📚 After Pipeline Completes

### 1. Explore Results
```python
import anndata as ad
import scanpy as sc

# Load
adata = ad.read_h5ad("results/analysis/final/bcell_atlas.h5ad")

# Inspect
print(f"Cells: {adata.n_obs:,}")
print(f"Genes: {adata.n_vars:,}")
print(f"Metadata: {list(adata.obs.columns)}")
print(f"Clusters: {adata.obs['leiden'].nunique()}")

# Visualize
sc.pl.umap(adata, color='leiden')
```

### 2. Quality Reports
```bash
# Open web reports
open results/analysis/final/cellranger_web_summary.html  # Cell Ranger QC
open results/report.html                                  # Nextflow execution report
open results/dag.svg                                      # Pipeline DAG
```

### 3. Cluster Annotation
```r
library(Seurat)
seurat <- as.Seurat(adata)

# Find markers
markers <- FindAllMarkers(seurat)

# Assign cell types based on markers
# (downstream analysis outside this pipeline)
```

## 🔗 Reference Links

- **Documentation**: [README.md](README.md) | [Quick Reference](QUICK_REFERENCE.md)
- **Step-by-step guides**: [STEP_1_COMPLETE.md](STEP_1_COMPLETE.md) | [STEP_2_COMPLETE.md](STEP_2_COMPLETE.md) | [STEP_3_COMPLETE.md](STEP_3_COMPLETE.md)
- **Container setup**: [CONTAINERS_GUIDE.md](CONTAINERS_GUIDE.md)
- **Nextflow**: https://www.nextflow.io/docs/latest/
- **10x Genomics**: https://support.10xgenomics.com/
- **Scanpy**: https://scanpy.readthedocs.io/
- **AnnData**: https://anndata.readthedocs.io/

## 🆘 Getting Help

1. **Check logs**:
   ```bash
   tail -f .nextflow.log
   tail -f results/trace.txt
   ```

2. **Review failed process**:
   ```bash
   # Check work directory for failed task
   ls results/work/
   cat results/work/*/command.log
   ```

3. **Consult troubleshooting**:
   - See [QUICK_REFERENCE.md](QUICK_REFERENCE.md#troubleshooting) for common issues
   - Check [CONTAINERS_GUIDE.md](CONTAINERS_GUIDE.md#troubleshooting) for container problems

4. **File an issue**:
   - Include error messages from logs
   - Provide sample sheet (without sensitive paths)
   - Describe system setup (CPU, GPU, OS)

---

## 🎯 Summary: You're Ready to Start!

**Next steps:**
1. ✅ Review this checklist
2. ✅ Install dependencies (Nextflow, Singularity)
3. ✅ Build containers: `bash scripts/build_containers.sh`
4. ✅ Prepare sample sheet (my_samples.csv)
5. ✅ Run test: `./nextflow run main.nf -profile desktop --sample_sheet test_data/sample_sheet.csv`
6. ✅ Process your data: `./nextflow run main.nf -profile desktop --sample_sheet my_samples.csv --outdir results`

**Expected runtime:**
- 10K cells: ~15 minutes
- 100K cells: ~1 hour
- 1M cells: ~5 hours

**Questions?** See [QUICK_REFERENCE.md](QUICK_REFERENCE.md) or check the troubleshooting guides.

Good luck! 🚀
