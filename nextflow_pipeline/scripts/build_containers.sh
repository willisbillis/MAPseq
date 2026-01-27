#!/bin/bash
# scripts/build_containers.sh
# Build and cache Singularity containers for the pipeline

set -e

# Configuration
CONTAINER_CACHE="${1:=${HOME}/.singularity/cache}"
CONTAINER_DEF_DIR="${2:=containers}"
QUIET="${3:-false}"

# Container specifications
declare -A CONTAINERS=(
    [rapids-cellranger]="rapids-cellranger.def"
    [analysis-gpu]="analysis-gpu.def"
)

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║          Singularity Container Build Script                   ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "Configuration:"
echo "  • Cache directory: $CONTAINER_CACHE"
echo "  • Definition files: $CONTAINER_DEF_DIR"
echo "  • Quiet mode: $QUIET"
echo ""

# Create cache directory
mkdir -p "$CONTAINER_CACHE"
echo "✓ Cache directory ready: $CONTAINER_CACHE"
echo ""

# Build containers
for container_name in "${!CONTAINERS[@]}"; do
    def_file="${CONTAINER_DEF_DIR}/${CONTAINERS[$container_name]}"
    output_sif="${CONTAINER_CACHE}/${container_name}-latest.sif"
    
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Building: $container_name"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    # Check definition file
    if [ ! -f "$def_file" ]; then
        echo "✗ ERROR: Definition file not found: $def_file"
        exit 1
    fi
    
    echo "Definition: $def_file"
    echo "Output: $output_sif"
    echo ""
    
    # Check if already exists
    if [ -f "$output_sif" ]; then
        echo "Container already exists. Remove to rebuild:"
        echo "  rm $output_sif"
        ls -lh "$output_sif"
        echo ""
        continue
    fi
    
    # Build container
    echo "Building container (this may take 10-30 minutes)..."
    echo ""
    
    if sudo singularity build "$output_sif" "$def_file"; then
        echo ""
        echo "✓ Build successful"
        ls -lh "$output_sif"
        
        # Test container
        echo ""
        echo "Testing container..."
        singularity exec "$output_sif" python3 -c "
import sys
print('Python version:', sys.version)
try:
    import scanpy as sc
    print('✓ Scanpy:', sc.__version__)
except ImportError:
    print('⚠ Scanpy not available')
"
    else
        echo ""
        echo "✗ Build failed for $container_name"
        exit 1
    fi
    
    echo ""
done

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║              Container Build Complete                         ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "Summary:"
echo "  • Cache directory: $CONTAINER_CACHE"
echo "  • Containers built: $(ls -1 $CONTAINER_CACHE/*.sif 2>/dev/null | wc -l)"
echo ""
echo "To use containers in Nextflow:"
echo "  export SINGULARITY_CACHE_DIR=\"$CONTAINER_CACHE\""
echo "  nextflow run main.nf -profile desktop"
echo ""
