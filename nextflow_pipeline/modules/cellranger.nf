// modules/cellranger.nf
// Cell Ranger multi orchestration with dynamic multi_config.csv generation

process CELLRANGER_MULTI_PREPARE {
    label 'cpu_small'
    publishDir "${params.outdir}/cellranger_configs", mode: 'copy'
    
    input:
    val(all_samples)  // List of sample maps from sample sheet
    
    output:
    path("multi_config_*.csv"), emit: configs
    path("config_generation_report.txt"), emit: report
    
    script:
    """
    python3 << 'PYSCRIPT'
    import csv
    import json
    import os
    from pathlib import Path
    from collections import defaultdict
    import sys
    
    print("=" * 70)
    print("Generating Cell Ranger multi_config.csv Files")
    print("=" * 70)
    
    # Parse samples data
    samples_data = ${all_samples}
    
    if not samples_data:
        print("ERROR: No samples provided")
        sys.exit(1)
    
    # Group samples by batch
    batches = defaultdict(list)
    for sample in samples_data:
        batches[sample['batch_id']].append(sample)
    
    print(f"\\nTotal batches: {len(batches)}")
    print(f"Total samples: {len(samples_data)}")
    
    # Track statistics
    stats = {
        'batches': len(batches),
        'samples_per_batch': {},
        'files_generated': 0,
        'vdj_count': 0,
        'hto_count': 0
    }
    
    # Generate multi_config.csv for each batch
    for batch_id, batch_samples in sorted(batches.items()):
        config_path = f"multi_config_{batch_id}.csv"
        
        print(f"\\nBatch: {batch_id}")
        print(f"  Samples: {len(batch_samples)}")
        
        with open(config_path, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # [gene-expression] section
            writer.writerow(['[gene-expression]'])
            writer.writerow(['reference', "${params.gex_ref}"])
            writer.writerow(['create-bam', 'true'])
            writer.writerow([])
            
            # [feature] section (HTO/ADT)
            has_hto = any(s['hto_r1'] is not None for s in batch_samples)
            if has_hto:
                writer.writerow(['[feature]'])
                writer.writerow(['reference', "${params.hto_feature_ref}"])
                writer.writerow([])
            
            # [vdj] section (if present)
            has_vdj = any(s['vdj_r1'] is not None for s in batch_samples)
            if has_vdj:
                writer.writerow(['[vdj]'])
                writer.writerow(['reference', "${params.vdj_ref}"])
                writer.writerow([])
            
            # [libraries] section
            writer.writerow(['[libraries]'])
            writer.writerow(['fastq_id,fastqs,feature_types'])
            
            lib_count = 0
            for sample in batch_samples:
                # GEX library (required)
                if sample['gex_r1'] and sample['gex_r2']:
                    fastq_dir = os.path.dirname(sample['gex_r1'])
                    writer.writerow([
                        f"{sample['sample_id']}_GEX",
                        fastq_dir,
                        "gene expression"
                    ])
                    lib_count += 1
                
                # HTO/ADT library
                if sample['hto_r1'] and sample['hto_r2']:
                    fastq_dir = os.path.dirname(sample['hto_r1'])
                    writer.writerow([
                        f"{sample['sample_id']}_HTO",
                        fastq_dir,
                        "antibody capture"
                    ])
                    lib_count += 1
                    stats['hto_count'] += 1
                
                # VDJ library
                if sample['vdj_r1'] and sample['vdj_r2']:
                    fastq_dir = os.path.dirname(sample['vdj_r1'])
                    vdj_type = sample.get('vdj_type', 'vdj-bcr')
                    writer.writerow([
                        f"{sample['sample_id']}_VDJ",
                        fastq_dir,
                        vdj_type
                    ])
                    lib_count += 1
                    stats['vdj_count'] += 1
        
        print(f"  Libraries: {lib_count}")
        print(f"  VDJ samples: {sum(1 for s in batch_samples if s['vdj_r1'])}")
        print(f"  HTO samples: {sum(1 for s in batch_samples if s['hto_r1'])}")
        print(f"  ✓ Generated: {config_path}")
        
        stats['samples_per_batch'][batch_id] = len(batch_samples)
        stats['files_generated'] += 1
    
    # Generate report
    with open("config_generation_report.txt", "w") as f:
        f.write("Cell Ranger multi_config.csv Generation Report\\n")
        f.write("=" * 70 + "\\n\\n")
        f.write(f"Batches: {stats['batches']}\\n")
        f.write(f"Total samples: {sum(stats['samples_per_batch'].values())}\\n")
        f.write(f"Files generated: {stats['files_generated']}\\n\\n")
        f.write(f"Modality summary:\\n")
        f.write(f"  • VDJ samples: {stats['vdj_count']}\\n")
        f.write(f"  • HTO samples: {stats['hto_count']}\\n\\n")
        f.write(f"Batches:\\n")
        for batch_id, count in stats['samples_per_batch'].items():
            f.write(f"  • {batch_id}: {count} samples\\n")
    
    print("\\n" + "=" * 70)
    print("✓ Config generation completed")
    print("=" * 70)
    PYSCRIPT
    """
}

process CELLRANGER_MULTI_EXECUTE {
    label 'cpu_large'
    publishDir "${params.outdir}/cellranger_outs/${batch_id}", mode: 'move', pattern: "*/outs/**"
    
    input:
    path(multi_config)  // multi_config.csv for a specific batch
    
    output:
    tuple val(batch_id), val(sample_ids), path("*/outs"), emit: outputs
    path("*/outs/web_summary.html"), optional: true, emit: web_summary
    path("cellranger_log_${batch_id}.txt"), emit: log
    
    script:
    batch_id = multi_config.baseName.replaceAll("multi_config_", "")
    sample_ids = "batch_${batch_id}"  // Placeholder, will be enriched by post-processing
    
    """
    #!/bin/bash
    set -e
    
    echo "╔════════════════════════════════════════════════════════════════╗"
    echo "║  Cell Ranger Multi: ${batch_id}                                  ║"
    echo "╚════════════════════════════════════════════════════════════════╝"
    
    LOG_FILE="cellranger_log_${batch_id}.txt"
    > "\$LOG_FILE"
    
    echo "Batch: ${batch_id}" | tee -a "\$LOG_FILE"
    echo "Config: ${multi_config}" | tee -a "\$LOG_FILE"
    echo "Start time: \$(date)" | tee -a "\$LOG_FILE"
    echo "" | tee -a "\$LOG_FILE"
    
    # Create working directory
    mkdir -p cellranger_${batch_id}
    cd cellranger_${batch_id}
    
    # Copy config
    cp ../${multi_config} multi_config.csv
    
    # Log config content for debugging
    echo "Config content:" >> "\$LOG_FILE"
    cat multi_config.csv | head -20 >> "\$LOG_FILE"
    echo "" >> "\$LOG_FILE"
    
    # Execute cellranger multi with error handling
    echo "Running cellranger multi..." | tee -a "\$LOG_FILE"
    
    if cellranger multi \\
        --id=run_${batch_id} \\
        --csv=multi_config.csv \\
        --localcores=${params.cellranger_localcores} \\
        --localmem=${params.cellranger_localmem} \\
        2>&1 | tee -a "\$LOG_FILE"; then
        
        echo "✓ Cell Ranger completed successfully" | tee -a "\$LOG_FILE"
        
        # Verify critical outputs
        if [ -d "run_${batch_id}/outs" ]; then
            OUTPUT_SIZE=\$(du -sh run_${batch_id}/outs | awk '{print \$1}')
            echo "✓ Output directory: run_${batch_id}/outs (\$OUTPUT_SIZE)" | tee -a "\$LOG_FILE"
            
            # List key outputs
            ls -lh run_${batch_id}/outs/ | grep -E "\\.(h5|csv|bam)$" | head -5 >> "\$LOG_FILE"
        else
            echo "✗ ERROR: No output directory found" | tee -a "\$LOG_FILE"
            exit 1
        fi
    else
        echo "✗ Cell Ranger failed with exit code \$?" | tee -a "\$LOG_FILE"
        exit 1
    fi
    
    # Copy web summary and log to parent
    cd ..
    cp cellranger_${batch_id}/run_${batch_id}/outs/web_summary.html run_${batch_id}_web_summary.html 2>/dev/null || true
    cp cellranger_${batch_id}/"\$LOG_FILE" .
    
    echo "End time: \$(date)" | tee -a "\$LOG_FILE"
    """
}

process CELLRANGER_COUNTS_TO_H5AD {
    label 'cpu_medium'
    publishDir "${params.outdir}/cellranger_h5ad/${batch_id}", mode: 'copy'
    
    input:
    tuple val(batch_id), path(cellranger_outs)
    
    output:
    tuple val(batch_id), path("raw_counts.h5ad"), emit: h5ad
    path("conversion_report.txt"), emit: report
    
    script:
    """
    python3 << 'PYSCRIPT'
    import scanpy as sc
    import numpy as np
    from pathlib import Path
    import sys
    
    print("Converting Cell Ranger output to AnnData H5AD")
    print("=" * 60)
    
    # Find Cell Ranger HDF5 output
    cr_outs = "${cellranger_outs}"
    h5_candidates = [
        Path(cr_outs) / "filtered_feature_bc_matrix.h5",
        Path(cr_outs) / "filtered_matrix.h5"
    ]
    
    cr_h5 = None
    for candidate in h5_candidates:
        if candidate.exists():
            cr_h5 = str(candidate)
            break
    
    if not cr_h5:
        print(f"ERROR: Cell Ranger output not found in {cr_outs}")
        print(f"Candidates: {h5_candidates}")
        sys.exit(1)
    
    print(f"Input: {cr_h5}")
    
    # Read count matrix
    try:
        adata = sc.read_h5ad(cr_h5)
    except Exception as e:
        print(f"ERROR reading H5AD: {e}")
        sys.exit(1)
    
    print(f"Cells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")
    
    # Add batch metadata
    adata.obs['batch'] = "${batch_id}"
    
    # Initial QC metrics
    adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()
    adata.obs['n_genes'] = np.array((adata.X > 0).sum(axis=1)).flatten()
    
    # Store Cell Ranger metadata
    adata.uns['cellranger'] = {
        'batch': "${batch_id}",
        'output_path': cr_outs
    }
    
    # Save to H5AD
    output_path = "raw_counts.h5ad"
    adata.write_h5ad(output_path, compression="gzip")
    print(f"✓ Saved: {output_path}")
    
    # Generate report
    with open("conversion_report.txt", "w") as f:
        f.write(f"Cell Ranger Conversion Report\\n")
        f.write(f"Batch: ${batch_id}\\n")
        f.write(f"Cells: {adata.n_obs}\\n")
        f.write(f"Genes: {adata.n_vars}\\n")
        f.write(f"Mean UMI per cell: {adata.obs['n_counts'].mean():.0f}\\n")
        f.write(f"Median UMI per cell: {adata.obs['n_counts'].median():.0f}\\n")
        f.write(f"Mean genes per cell: {adata.obs['n_genes'].mean():.0f}\\n")
        f.write(f"Output: {output_path}\\n")
    PYSCRIPT
    """
}
