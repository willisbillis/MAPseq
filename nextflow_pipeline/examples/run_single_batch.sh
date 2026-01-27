#!/bin/bash
# examples/run_single_batch.sh
# Minimal example: single batch on desktop

# Default values (override via command line)
SAMPLE_SHEET="${1:-test_data/sample_sheet.csv}"
REF_DIR="${2:-/home/Apps/genomes/cellranger}"
OUTDIR="${3:-./results_single_batch}"
GENOME="${4:-GRCh38}"

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║  B-Cell Atlas Pipeline - Single Batch Example (Desktop)       ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "Configuration:"
echo "  • Sample sheet: $SAMPLE_SHEET"
echo "  • Reference dir: $REF_DIR"
echo "  • Output dir: $OUTDIR"
echo "  • Genome: $GENOME"
echo ""

# Run Nextflow pipeline
nextflow run main.nf \
    -profile desktop \
    -resume \
    --sample_sheet "$SAMPLE_SHEET" \
    --ref_dir "$REF_DIR" \
    --genome "$GENOME" \
    --outdir "$OUTDIR" \
    --max_parallel_batches 1 \
    --cellbender_enabled true \
    --harmony_enabled true \
    --include_vdj true \
    --hca_compliance_enabled true

echo ""
echo "✓ Pipeline completed!"
echo "Results saved to: $OUTDIR"
