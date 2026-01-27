#!/usr/bin/env nextflow
// main.nf
// Master workflow orchestration for B-cell atlas pipeline

nextflow.enable.dsl = 2

// Workflow title
println """
    ╔══════════════════════════════════════════════════════════════════════════════╗
    ║          B-Cell Atlas Nextflow Pipeline v${workflow.manifest.version}                    ║
    ║                GPU-accelerated single-cell analysis                          ║
    ╚══════════════════════════════════════════════════════════════════════════════╝
    
    Inputs:
      • Sample sheet:    ${params.sample_sheet}
      • Reference dir:   ${params.ref_dir}
      • Output dir:      ${params.outdir}
      • Genome:          ${params.genome}
      • GPU strategy:    ${params.gpu_strategy}
    
    Processing:
      • CellBender QC:   ${params.cellbender_enabled ? 'ENABLED' : 'DISABLED'}
      • Harmony batch:   ${params.harmony_enabled ? 'ENABLED' : 'DISABLED'}
      • VDJ analysis:    ${params.include_vdj == null ? 'AUTO-DETECT' : params.include_vdj ? 'ENABLED' : 'DISABLED'}
      • HCA compliance:  ${params.hca_compliance_enabled ? 'ENABLED' : 'DISABLED'}
    
    Resources:
      • Workflow executor: ${workflow.executor}
      • Available threads: ${Runtime.getRuntime().availableProcessors()}
      • Max memory: ${Runtime.getRuntime().maxMemory() / 1024 / 1024 / 1024} GB
      • Start time: ${workflow.start}
"""

// Include workflow modules
include { PREFLIGHT_CHECKS } from './modules/preflight.nf'
include { CELLRANGER_MULTI_PREPARE } from './modules/cellranger.nf'
include { CELLRANGER_MULTI_EXECUTE } from './modules/cellranger.nf'
include { AGGREGATION_PREPARE } from './modules/aggregation.nf'
include { AGGREGATION_EXECUTE } from './modules/aggregation.nf'
include { CELLBENDER_QC; SCDBLFINDER_FILTER } from './modules/qc.nf'
include { HARMONY_INTEGRATION } from './modules/analysis.nf'
include { RAPIDS_CLUSTERING } from './modules/analysis.nf'
include { METADATA_INJECT; METADATA_VALIDATE } from './modules/metadata.nf'
include { VDJ_CLONOTYPE_ANALYSIS; VDJ_INTEGRATION } from './modules/vdj.nf'
include { PORTAL_MANIFEST } from './modules/portal_manifest.nf'

// ============================================================================
// CHANNEL SETUP & DATA LOADING
// ============================================================================

workflow {
    
    // Step 0: Preflight checks
    PREFLIGHT_CHECKS()
    
    // Step 1: Load sample sheet and parse into channels
    samples_ch = Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row ->
            [
                sample_id: row.sample_id,
                batch_id: row.batch_id,
                gex_r1: row.gex_fastq_r1 ?: null,
                gex_r2: row.gex_fastq_r2 ?: null,
                hto_r1: row.hto_fastq_r1 ?: null,
                hto_r2: row.hto_fastq_r2 ?: null,
                vdj_r1: row.vdj_fastq_r1 ?: null,
                vdj_r2: row.vdj_fastq_r2 ?: null,
                donor_id: row.donor_id ?: 'unknown',
                tissue: row.tissue ?: params.organ,
                facs_gates: row.facs_gates ?: null
            ]
        }
    
    // Auto-detect VDJ presence in sample sheet if not specified
    vdj_present = samples_ch
        .map { it.vdj_r1 != null && it.vdj_r2 != null }
        .collect()
        .map { list -> list.any() }
    
    include_vdj_flag = params.include_vdj != null ? params.include_vdj : vdj_present.first()
    
    if (params.debug) {
        println "DEBUG: VDJ auto-detect = ${vdj_present.first()}, include_vdj flag = ${include_vdj_flag}"
    }
    
    // Step 2: Cellranger multi_config.csv generation
    cellranger_multi_configs = CELLRANGER_MULTI_PREPARE(
        samples_ch.collect()
    )
    
    // Step 3: Execute cellranger multi per batch (parallel)
    cellranger_outputs = CELLRANGER_MULTI_EXECUTE(
        cellranger_multi_configs.flatten()
    )
    
    // Step 4: Prepare aggregation manifests per batch
    aggregation_manifests = AGGREGATION_PREPARE(
        cellranger_outputs.groupTuple(by: 1)  // Group by batch_id
    )
    
    // Step 5: Execute aggregation per batch
    aggregated_matrices = AGGREGATION_EXECUTE(
        aggregation_manifests
    )
    
    // Step 6: Quality control (CellBender + scDblFinder)
    qc_filtered = CELLBENDER_QC(
        aggregated_matrices,
        PREFLIGHT_CHECKS.out.gpu_available
    )
    
    qc_doublets = SCDBLFINDER_FILTER(
        qc_filtered
    )
    
    // Step 7: Metadata validation and injection
    metadata_validated = METADATA_VALIDATE(
        qc_doublets
    )
    
    metadata_injected = METADATA_INJECT(
        metadata_validated,
        samples_ch.collect()
    )
    
    // Step 8: Batch integration (Harmony)
    if (params.harmony_enabled) {
        harmony_output = HARMONY_INTEGRATION(
            metadata_injected,
            PREFLIGHT_CHECKS.out.gpu_available
        )
        pre_clustering = harmony_output
    } else {
        pre_clustering = metadata_injected
    }
    
    // Step 9: Clustering and dimensional reduction (rapids-singlecell)
    clustering_output = RAPIDS_CLUSTERING(
        pre_clustering,
        PREFLIGHT_CHECKS.out.gpu_available
    )
    
    // Step 10: Optional VDJ clonotype analysis (conditional)
    if (include_vdj_flag) {
        vdj_stats = VDJ_CLONOTYPE_ANALYSIS(
            cellranger_outputs
                .groupTuple(by: 1)
                .map { it[2] }  // Extract VDJ outputs
        )
        
        final_output = VDJ_INTEGRATION(
            clustering_output,
            vdj_stats
        )
    } else {
        final_output = clustering_output
        println "INFO: Skipping VDJ analysis (no VDJ data detected or --include_vdj=false)"
    }

    // Step 11b: Generate portal-ready manifest (JSON + CSV)
    PORTAL_MANIFEST(
        final_output
    )
    
    // Step 11: Emit final outputs
    final_output
}

// ============================================================================
// WORKFLOW COMPLETION REPORT
// ============================================================================

workflow.onComplete {
    println """
    ╔══════════════════════════════════════════════════════════════════════════════╗
    ║                      Pipeline Execution Completed                            ║
    ╚══════════════════════════════════════════════════════════════════════════════╝
    
    Duration:       ${workflow.complete ? workflow.duration : 'N/A'}
    Success:        ${workflow.success ? 'YES ✓' : 'FAILED ✗'}
    Exit status:    ${workflow.exitStatus}
    
    Results saved to: ${params.outdir}
    
    Key outputs:
      • Primary analysis:  ${params.outdir}/analysis/bcell_atlas.h5ad
      • Checkpoints:       ${params.outdir}/checkpoints/
      • Reports:           ${params.outdir}/trace.txt, timeline.html
    
    ${workflow.success ? '✓ Pipeline completed successfully!' : '✗ Pipeline failed. Check logs in ' + params.outdir}
    """
}

workflow.onError {
    println """
    ╔══════════════════════════════════════════════════════════════════════════════╗
    ║                        Pipeline Execution Failed                             ║
    ╚══════════════════════════════════════════════════════════════════════════════╝
    
    Error message:
      ${workflow.errorMessage ?: 'Unknown error'}
    
    Check the following for debugging:
      • Log files: ${params.outdir}/
      • Trace report: ${params.outdir}/trace.txt
      • DAG visualization: ${params.outdir}/dag.svg
    """
}
