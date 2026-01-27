# Singularity Container Setup Guide

## Overview

The B-cell atlas Nextflow pipeline uses two specialized Singularity containers:

1. **rapids-cellranger.def** — Cell Ranger orchestration + GPU support
   - Base: NVIDIA RAPIDS 25.12 (CUDA 12.0)
   - Tools: Cell Ranger 8.0, Python analysis libs
   - Use case: cellranger multi execution

2. **analysis-gpu.def** — GPU-accelerated analysis
   - Base: NVIDIA RAPIDS 25.12 (CUDA 12.0)
   - Tools: Scanpy, rapids-singlecell, Harmony (PyTorch), scirpy
   - Use case: QC, clustering, batch integration, VDJ analysis

## Quick Start

### 1. Install Singularity

```bash
# Ubuntu/Debian
sudo apt-get install -y singularity-container

# macOS (requires Docker Desktop or Multipass)
brew install --cask multipass singularity

# HPC clusters: typically already installed
module load singularity
singularity --version
```

### 2. Build Containers

```bash
# Desktop/Linux
bash scripts/build_containers.sh ~/.singularity/cache containers/

# With custom cache directory (HPC)
bash scripts/build_containers.sh /mnt/shared/singularity/cache containers/

# Build individual container
sudo singularity build ~/.singularity/cache/rapids-cellranger-latest.sif containers/rapids-cellranger.def
```

**Build time:** ~15-30 minutes per container (depending on network speed and hardware)

**Output:**
```
~/.singularity/cache/
├── rapids-cellranger-latest.sif    (~3 GB)
└── analysis-gpu-latest.sif         (~2.5 GB)
```

### 3. Verify Installation

```bash
# Test Cell Ranger container
singularity exec --nv ~/.singularity/cache/rapids-cellranger-latest.sif python3 -c "import scanpy; print('✓ Scanpy installed')"

# Test analysis container with GPU
singularity exec --nv ~/.singularity/cache/analysis-gpu-latest.sif python3 << 'EOF'
import torch
print(f"GPU available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"GPU: {torch.cuda.get_device_name(0)}")
EOF
```

## Container Details

### rapids-cellranger.def

**Base image:** `nvcr.io/nvidia/rapids:25.12-cuda12.0-runtime-ubuntu22.04`

**Key packages:**
- **Cell Ranger 8.0.0** (placeholder; must be installed separately or via binary)
- **Scanpy 1.10.0**
- **rapids-singlecell** (GPU acceleration)
- **cuML 25.12** (GPU machine learning)
- **harmonypy 0.2.0+** (PyTorch GPU batch integration)
- **scirpy 0.13** (VDJ/TCR analysis)
- **R 4.2** (optional Seurat support)

**GPU support:**
- NVIDIA CUDA 12.0
- cuDNN integrated
- PyTorch with CUDA backend

**Usage in Nextflow:**
```groovy
container = "${params.singularity_cache_dir}/rapids-cellranger-latest.sif"
containerOptions = '--nv --gpus=0'  // Bind GPU device 0
```

### analysis-gpu.def

**Base image:** `nvcr.io/nvidia/rapids:25.12-cuda12.0-runtime-ubuntu22.04`

**Key packages:**
- **Scanpy 1.10.0** (single-cell analysis)
- **rapids-singlecell 0.12.0+** (GPU PCA, UMAP, Leiden)
- **harmonypy 0.2.0+** (GPU batch integration with PyTorch backend)
- **scirpy 0.13.0** (V(D)J repertoire analysis)
- **PyTorch 2.0+** (CUDA backend for Harmony)
- **Optional:** CellBender, SoupX, doubletdetection

**GPU support:**
- CUDA 12.0 runtime
- PyTorch 2.0+ with CUDA
- rapids-singlecell GPU acceleration
- cuML GPU algorithms

**Usage in Nextflow:**
```groovy
container = "${params.singularity_cache_dir}/analysis-gpu-latest.sif"
containerOptions = '--nv --gpus=0'  // GPU-enabled analysis
```

## Advanced Configuration

### Using Containers on HPC

#### SLURM with GPU allocation

```bash
# In conf/hpc.config
singularity {
    cacheDir = "/mnt/shared/singularity/cache"
    runOptions = "--bind /mnt/shared,/scratch --nv"
}

process {
    withLabel: 'gpu_.*' {
        clusterOptions = "--gres=gpu:1"
        containerOptions = '--nv --gpus=${task.gpu}'
    }
}
```

**Build on HPC head node:**
```bash
# Request interactive GPU node
salloc --gpus=1 --time=1:00:00 -n 4 --mem=16G

# Build container
singularity build /mnt/shared/singularity/cache/analysis-gpu-latest.sif containers/analysis-gpu.def
```

#### Pre-cached containers

For clusters with many nodes, pre-cache containers on shared storage:

```bash
# Head node (once)
bash scripts/build_containers.sh /mnt/shared/singularity/cache containers/

# All compute nodes: containers automatically mounted via NFS
# Verify:
ls -la /mnt/shared/singularity/cache/*.sif
```

### Using Custom Cell Ranger Binary

The `rapids-cellranger.def` container does not include the Cell Ranger binary (it's large and requires 10x Genomics licensing). To include Cell Ranger:

**Option 1: Mount host Cell Ranger**
```bash
# If cellranger installed on host at /opt/cellranger
nextflow run main.nf --singularity_container_options "--bind /opt/cellranger:/cellranger"

# In params.config:
cellranger_cmd = "singularity exec --bind /opt/cellranger:/cellranger --nv ${params.container_cellranger} /cellranger/cellranger"
```

**Option 2: Custom container with Cell Ranger embedded**

Create `containers/rapids-cellranger-custom.def`:
```
Bootstrap: docker
From: nvcr.io/nvidia/rapids:25.12-cuda12.0-runtime-ubuntu22.04

%post
    # ... (install packages from rapids-cellranger.def)
    
    # Download and install Cell Ranger
    cd /opt
    wget https://cf.10xgenomics.com/supp/cell-exp/cellranger-8.0.0.tar.gz
    tar -xzf cellranger-8.0.0.tar.gz
    export PATH=/opt/cellranger-8.0.0:$PATH
```

### GPU Device Management

#### Auto-detection (recommended)

```bash
# Nextflow automatically detects available GPUs
nextflow run main.nf \
    --gpu_strategy auto \
    -profile slurm

# Automatically scales GPU usage:
# • 1-100K cells: 1 GPU (6GB VRAM)
# • 100K-500K cells: 2 GPUs (12GB VRAM)
# • 500K-1M cells: 4 GPUs (24GB VRAM)
```

#### Manual GPU assignment

```bash
# Force specific GPU
nextflow run main.nf \
    --gpu_devices_available 0,1 \
    -profile local

# In process:
containerOptions = '--nv --gpus=1,2'  # Use GPUs 1 and 2
```

## Troubleshooting

### Build fails: "NVIDIA base image not accessible"

**Solution:**
```bash
# Login to NVIDIA container registry
singularity remote login --username $oauthtoken docker://nvcr.io

# Or set credentials file
export SINGULARITY_DOCKER_LOGIN_FILE=~/.singularity/docker_login.json
```

### GPU not detected in container

```bash
# Verify host GPU
nvidia-smi

# Verify Singularity GPU support
singularity exec --nv --gpu /usr/bin/nvidia-smi

# If not working, update Singularity:
singularity --version  # Should be 3.7+
sudo apt-get upgrade singularity
```

### Container too large for cache

```bash
# Check cache size
du -sh ~/.singularity/cache/

# Use larger cache location
export SINGULARITY_CACHEDIR=/mnt/large_storage/.singularity/cache
mkdir -p $SINGULARITY_CACHEDIR
bash scripts/build_containers.sh $SINGULARITY_CACHEDIR containers/
```

### Python package conflicts

```bash
# Rebuild from scratch (clears cached layers)
rm ~/.singularity/cache/*.sif
sudo singularity build --force ~/.singularity/cache/analysis-gpu-latest.sif containers/analysis-gpu.def

# Or skip build cache:
sudo singularity build --no-cache ~/.singularity/cache/analysis-gpu-latest.sif containers/analysis-gpu.def
```

## Production Checklist

- ✅ Containers built and cached locally/on shared storage
- ✅ NVIDIA GPUs verified available (`nvidia-smi` works in container)
- ✅ CUDA 12.0 runtime compatible with host NVIDIA driver (≥520)
- ✅ Disk space available for containers (~5-6 GB total)
- ✅ Singularity version ≥ 3.7 installed
- ✅ Network access to NVIDIA container registry (for base images)
- ✅ PyTorch GPU backend working (`torch.cuda.is_available()` = True)
- ✅ rapids-singlecell GPU acceleration verified
- ✅ Cell Ranger binary accessible (via mount or local install)

## Performance Notes

### GPU Memory Usage (rapids-singlecell)

| Operation | Dataset Size | Memory | Runtime |
|-----------|-------------|--------|---------|
| PCA (50 PCs) | 100K cells × 2K genes | ~2 GB | 10s |
| UMAP | 100K cells × 50 PCs | ~1 GB | 5s |
| Leiden clustering | 100K cells × 50 PCs | ~0.5 GB | 3s |
| Harmony integration | 1M cells × 50 PCs | ~8 GB | 45s |

**CPU equivalents (Scanpy):**
- PCA: ~2 minutes (50x slower)
- UMAP: ~1 minute (10x slower)
- Harmony: ~10 minutes (10x slower)

### Container Size Optimization

If storage is limited:

```bash
# Minimal container (Scanpy CPU-only, no R)
# Edit analysis-gpu.def:
# • Remove R installation
# • Remove optional tools (CellBender, SoupX)
# • Remove optional R packages
# Result: ~1.5 GB

# Skip rapids-singlecell (CPU-only)
# Remove rapids-singlecell installation
# Result: ~2 GB, CPU analysis only
```

## Container Registry Integration (Optional)

For team sharing, push containers to registry:

```bash
# Convert to OCI format for registry
singularity build --sandbox sandbox/ analysis-gpu-latest.sif
singularity remote add mylab oras://mylab.azurecr.io

# Push to Azure Container Registry
singularity push sandbox/ oras://mylab.azurecr.io/bcell-atlas:v1.0

# Pull on another machine
singularity pull oras://mylab.azurecr.io/bcell-atlas:v1.0
```

---

**Next steps:** After containers are built, proceed to Step 4 (QC pipeline) to enable actual data processing.
