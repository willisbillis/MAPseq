// modules/preflight.nf
// Validation and preflight checks for pipeline execution

process PREFLIGHT_CHECKS {
    label 'cpu_small'
    publishDir "${params.outdir}/logs", mode: 'copy'
    
    output:
    path("preflight_report.txt"), emit: report
    env GPU_AVAILABLE, emit: gpu_available
    
    script:
    """
    #!/bin/bash
    set -e
    
    echo "╔════════════════════════════════════════════════════════════════╗"
    echo "║           PREFLIGHT VALIDATION CHECKS                         ║"
    echo "╚════════════════════════════════════════════════════════════════╝"
    
    # Initialize report file
    REPORT="preflight_report.txt"
    > "\$REPORT"
    
    echo "Pipeline Execution Environment Report" >> "\$REPORT"
    echo "Date: \$(date)" >> "\$REPORT"
    echo "====================================" >> "\$REPORT"
    echo "" >> "\$REPORT"
    
    # 1. Check Nextflow version
    echo "Nextflow:" >> "\$REPORT"
    if command -v nextflow &> /dev/null; then
        NF_VERSION=\$(nextflow -version 2>&1 | grep "Nextflow" | awk '{print \$NF}')
        echo "✓ Nextflow version: \$NF_VERSION" | tee -a "\$REPORT"
    else
        echo "✗ Nextflow not found" | tee -a "\$REPORT"
    fi
    echo "" >> "\$REPORT"
    
    # 2. Check reference directories
    echo "Reference Paths:" >> "\$REPORT"
    echo "---------------" >> "\$REPORT"
    
    REFS_OK=true
    for ref_path in "${params.gex_ref}" "${params.vdj_ref}" "${params.feature_ref_dir}"; do
        if [ -d "\$ref_path" ]; then
            echo "✓ \$ref_path" | tee -a "\$REPORT"
        else
            echo "⚠ MISSING: \$ref_path (may be on compute node)" | tee -a "\$REPORT"
            REFS_OK=false
        fi
    done
    echo "" >> "\$REPORT"
    
    # 3. Check required tools
    echo "Required Tools:" >> "\$REPORT"
    echo "---------------" >> "\$REPORT"
    
    TOOLS=("cellranger" "python3" "singularity")
    TOOLS_OK=true
    for tool in "\${TOOLS[@]}"; do
        if command -v \$tool &> /dev/null; then
            VERSION=\$(\$tool --version 2>&1 | head -1 | cut -d' ' -f1-3)
            echo "✓ \$tool: \$VERSION" | tee -a "\$REPORT"
        else
            echo "✗ MISSING: \$tool" | tee -a "\$REPORT"
            TOOLS_OK=false
        fi
    done
    echo "" >> "\$REPORT"
    
    # 4. Check GPU availability
    echo "GPU Configuration:" >> "\$REPORT"
    echo "------------------" >> "\$REPORT"
    
    GPU_AVAILABLE=false
    GPU_COUNT=0
    
    if command -v nvidia-smi &> /dev/null; then
        GPU_COUNT=\$(nvidia-smi --list-gpus 2>/dev/null | wc -l)
        echo "✓ NVIDIA GPUs available: \$GPU_COUNT" | tee -a "\$REPORT"
        nvidia-smi --query-gpu=index,name,driver_version --format=csv,noheader 2>/dev/null | sed 's/^/  /' >> "\$REPORT"
        GPU_AVAILABLE=true
    else
        echo "⚠ No NVIDIA GPUs detected (GPU acceleration will be unavailable)" | tee -a "\$REPORT"
        GPU_AVAILABLE=false
    fi
    echo "" >> "\$REPORT"
    
    # 5. Check container runtime
    echo "Container Runtime:" >> "\$REPORT"
    echo "------------------" >> "\$REPORT"
    
    if command -v singularity &> /dev/null; then
        SIF_VERSION=\$(singularity --version 2>/dev/null)
        echo "✓ Singularity: \$SIF_VERSION" | tee -a "\$REPORT"
    else
        echo "✗ Singularity not found" | tee -a "\$REPORT"
    fi
    
    if [ -d "${params.singularity_cache_dir}" ]; then
        CACHE_SIZE=\$(du -sh "${params.singularity_cache_dir}" 2>/dev/null | awk '{print \$1}')
        echo "✓ Singularity cache: ${params.singularity_cache_dir} (\$CACHE_SIZE)" | tee -a "\$REPORT"
    else
        echo "⚠ Creating Singularity cache: ${params.singularity_cache_dir}" | tee -a "\$REPORT"
        mkdir -p "${params.singularity_cache_dir}"
    fi
    echo "" >> "\$REPORT"
    
    # 6. Check sample sheet format
    echo "Sample Sheet:" >> "\$REPORT"
    echo "-------------" >> "\$REPORT"
    
    if [ -f "${params.sample_sheet}" ]; then
        SAMPLE_COUNT=\$(tail -n +2 "${params.sample_sheet}" 2>/dev/null | wc -l)
        BATCHES=\$(tail -n +2 "${params.sample_sheet}" 2>/dev/null | cut -d',' -f2 | sort -u | wc -l)
        echo "✓ Samples found: \$SAMPLE_COUNT" | tee -a "\$REPORT"
        echo "✓ Batches: \$BATCHES" | tee -a "\$REPORT"
        
        # Check for required columns
        HEADER=\$(head -n 1 "${params.sample_sheet}")
        for col in sample_id batch_id gex_fastq_r1 gex_fastq_r2; do
            if [[ "\$HEADER" == *"\$col"* ]]; then
                echo "  ✓ \$col" >> "\$REPORT"
            else
                echo "  ✗ Column missing: \$col" | tee -a "\$REPORT"
            fi
        done
    else
        echo "✗ Sample sheet not found: ${params.sample_sheet}" | tee -a "\$REPORT"
    fi
    echo "" >> "\$REPORT"
    
    # 7. Check disk space
    echo "Disk Space:" >> "\$REPORT"
    echo "-----------" >> "\$REPORT"
    
    if [ -d "${params.outdir}" ] || [ -d "\$(dirname "${params.outdir}")" ]; then
        OUTPUT_DISK=\$(df -BG "${params.outdir}" 2>/dev/null | tail -1 | awk '{print \$4}' | sed 's/G$//')
        echo "✓ Available in output directory: \${OUTPUT_DISK:-Unknown} GB" | tee -a "\$REPORT"
    else
        echo "⚠ Output directory does not exist, creating: ${params.outdir}" | tee -a "\$REPORT"
        mkdir -p "${params.outdir}"
    fi
    echo "" >> "\$REPORT"
    
    # 8. Memory availability
    echo "System Memory:" >> "\$REPORT"
    echo "--------------" >> "\$REPORT"
    
    TOTAL_MEM=\$(free -h 2>/dev/null | awk 'NR==2 {print \$2}')
    AVAIL_MEM=\$(free -h 2>/dev/null | awk 'NR==2 {print \$7}')
    echo "✓ Total: \$TOTAL_MEM" | tee -a "\$REPORT"
    echo "✓ Available: \$AVAIL_MEM" | tee -a "\$REPORT"
    echo "" >> "\$REPORT"
    
    # Summary
    echo "Preflight Checks Completed: \$(date)" >> "\$REPORT"
    
    # Check critical failures
    if [ "\$TOOLS_OK" = false ] && [ "\$REFS_OK" = false ]; then
        echo "" | tee -a "\$REPORT"
        echo "⚠ WARNING: Some critical checks failed. Please verify configuration." | tee -a "\$REPORT"
    fi
    
    echo ""
    echo "✓ Preflight checks completed. Details: \$REPORT"
    echo ""
    
    # Export GPU status for downstream processes
    echo "\$GPU_AVAILABLE"
    """
}
