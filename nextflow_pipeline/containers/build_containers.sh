#!/bin/bash

set -e

echo "======================================"
echo "MAPseq Singularity Container Builder"
echo "======================================"

# Check for Singularity/Apptainer
if command -v singularity &> /dev/null; then
    SINGULARITY_CMD="singularity"
    echo "✓ Found Singularity: $(singularity version)"
elif command -v apptainer &> /dev/null; then
    SINGULARITY_CMD="apptainer"
    echo "✓ Found Apptainer: $(apptainer version)"
else
    echo "✗ ERROR: Neither Singularity nor Apptainer found."
    echo "  Install from: https://singularity-ce.readthedocs.io/ or https://apptainer.org/"
    exit 1
fi

# Verify CUDA is available (optional but recommended for GPU containers)
if command -v nvidia-smi &> /dev/null; then
    echo "✓ NVIDIA GPU detected:"
    nvidia-smi --query-gpu=name --format=csv,noheader | head -1
else
    echo "⚠ WARNING: No NVIDIA GPU detected. Containers will still build but GPU execution may fail."
fi

# Build containers
CONTAINERS=("rapids-cellranger.def" "analysis-gpu.def")
OUTPUT_DIR="./build"

mkdir -p "$OUTPUT_DIR"

for def_file in "${CONTAINERS[@]}"; do
    if [ ! -f "$def_file" ]; then
        echo "✗ ERROR: $def_file not found"
        exit 1
    fi
    
    container_name="${def_file%.def}"
    output_file="$OUTPUT_DIR/${container_name}.sif"
    
    echo ""
    echo "Building: $container_name"
    echo "  Definition: $def_file"
    echo "  Output: $output_file"
    echo "  This may take 10-20 minutes..."
    echo ""
    
    # Build with sudo if needed (Singularity often requires root)
    if [ "$EUID" -eq 0 ]; then
        $SINGULARITY_CMD build "$output_file" "$def_file"
    else
        echo "NOTE: Building without sudo. If build fails due to permissions, run:"
        echo "  sudo $SINGULARITY_CMD build $output_file $def_file"
        $SINGULARITY_CMD build "$output_file" "$def_file" || {
            echo "Build failed. Trying with sudo..."
            sudo $SINGULARITY_CMD build "$output_file" "$def_file"
        }
    fi
    
    if [ -f "$output_file" ]; then
        size=$(du -h "$output_file" | cut -f1)
        echo "✓ Successfully built: $output_file ($size)"
    else
        echo "✗ Build failed for $container_name"
        exit 1
    fi
done

echo ""
echo "======================================"
echo "Build Complete!"
echo "======================================"
echo ""
echo "Containers ready in: $OUTPUT_DIR"
echo ""
echo "Test a container with:"
echo "  $SINGULARITY_CMD exec --nv build/rapids-cellranger.sif cellranger version"
echo "  $SINGULARITY_CMD exec --nv build/analysis-gpu.sif python3 -c 'import rapids; print(rapids.__version__)'"
echo ""
