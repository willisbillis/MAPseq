#!/bin/bash
# examples/run_multi_batch.sh
# Multi-batch example with Harmony integration and VDJ on HPC

SAMPLE_SHEET="${1:-test_data/sample_sheet.csv}"
REF_DIR="${2:-/home/Apps/genomes/cellranger}"
OUTDIR="${3:-./results_multi_batch}"
GENOME="${4:-GRCh38}"

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║  B-Cell Atlas Pipeline - Multi-Batch Example (HPC)            ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "Configuration:"
echo "  • Sample sheet: $SAMPLE_SHEET"
echo "  • Reference dir: $REF_DIR"
echo "  • Output dir: $OUTDIR"
echo "  • Genome: $GENOME"
echo ""
echo "Advanced parameters:"
echo "  • Harmony resolution: 50 PCs"
echo "  • Leiden clustering: 1.2 resolution"
echo "  • HVGs: 3000"
echo "  • GPU acceleration: ENABLED"
echo ""

# Run Nextflow pipeline on HPC
nextflow run main.nf \
    -profile slurm \
    -resume \
    -with-trace \
    --sample_sheet "$SAMPLE_SHEET" \
    --ref_dir "$REF_DIR" \
    --genome "$GENOME" \
    --outdir "$OUTDIR" \
    --max_parallel_batches 4 \
    --cellbender_enabled true \
    --cellbender_epochs_max 100 \
    --harmony_enabled true \
    --harmony_npcs 50 \
    --n_hvg 3000 \
    --leiden_resolution 1.2 \
    --include_vdj true \
    --hca_compliance_enabled true \
    --gpu_strategy auto

echo ""
echo "✓ Pipeline queued on HPC!"
echo "Results will be saved to: $OUTDIR"
echo "Monitor with: tail -f $OUTDIR/trace.txt"
