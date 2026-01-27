// modules/aggregation.nf
// Per-batch aggregation of Cell Ranger outputs and preparation for QC

process AGGREGATION_PREPARE {
    label 'cpu_small'
    publishDir "${params.outdir}/aggregation_manifests/${batch_id}", mode: 'copy'
    
    input:
    tuple val(batch_id), val(sample_ids), path(cellranger_outs)
    
    output:
    tuple val(batch_id), path("aggregation.csv"), path(cellranger_outs), emit: manifest
    
    script:
    """
    #!/bin/bash
    
    # Create aggregation manifest for cellranger aggr command
    # Format: library_id,molecule_h5
    
    > aggregation.csv
    echo "library_id,molecule_h5" >> aggregation.csv
    
    for cr_out in ${cellranger_outs}; do
        if [ -f "\$cr_out/molecule_info.h5" ]; then
            LIBRARY_ID=\$(basename \$(dirname \$cr_out))
            echo "\$LIBRARY_ID,\$cr_out/molecule_info.h5" >> aggregation.csv
        fi
    done
    
    echo "Generated aggregation manifest:"
    cat aggregation.csv
    """
}

process AGGREGATION_EXECUTE {
    label 'cpu_large'
    publishDir "${params.outdir}/aggregated/${batch_id}", mode: 'move', pattern: "*/outs/**"
    
    input:
    tuple val(batch_id), path(manifest), path(cellranger_outs)
    
    output:
    tuple val(batch_id), path("aggr_${batch_id}/outs"), emit: aggregated
    path("aggr_${batch_id}/outs/web_summary.html"), optional: true, emit: web_summary
    
    script:
    """
    #!/bin/bash
    set -e
    
    echo "================================================"
    echo "CellRanger Aggregation: ${batch_id}"
    echo "================================================"
    
    # Execute cellranger aggr
    cellranger aggr \\
        --id=aggr_${batch_id} \\
        --csv=${manifest} \\
        --normalize=mapped \\
        --localcores=${params.cellranger_localcores} \\
        --localmem=${params.cellranger_localmem}
    
    # Verify output
    if [ -d "aggr_${batch_id}/outs" ]; then
        echo "✓ Aggregation completed successfully"
        
        # List aggregated outputs
        echo "Generated outputs:"
        ls -lh aggr_${batch_id}/outs/ | grep -E "\\.(h5|csv)$"
    else
        echo "✗ Aggregation failed"
        exit 1
    fi
    """
}

process AGGREGATION_TO_H5AD {
    label 'cpu_medium'
    publishDir "${params.outdir}/h5ad_aggregated/${batch_id}", mode: 'copy'
    
    input:
    tuple val(batch_id), path(aggr_outs)
    
    output:
    tuple val(batch_id), path("aggregated_raw.h5ad"), emit: h5ad
    path("aggregation_stats.txt"), emit: stats
    
    script:
    """
    python3 << 'PYSCRIPT'
    import scanpy as sc
    import numpy as np
    from pathlib import Path
    import json
    
    # Load aggregated count matrix
    h5_file = "${aggr_outs}/filtered_feature_bc_matrix.h5"
    
    if not Path(h5_file).exists():
        raise FileNotFoundError(f"Aggregated matrix not found: {h5_file}")
    
    print(f"Loading aggregated count matrix from: {h5_file}")
    adata = sc.read_h5ad(h5_file)
    
    # Add batch information from metadata
    # This will be enriched in the metadata injection step
    adata.obs['batch'] = "${batch_id}"
    
    # Initial QC metrics (pre-filtering)
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1 if hasattr(adata.X, 'A1') else adata.X.sum(axis=1)
    adata.obs['n_genes'] = (adata.X > 0).sum(axis=1).A1 if hasattr(adata.X, 'A1') else (adata.X > 0).sum(axis=1)
    
    # Save aggregated object
    adata.write_h5ad("aggregated_raw.h5ad", compression="gzip")
    
    # Generate statistics
    with open("aggregation_stats.txt", "w") as f:
        f.write(f"Batch ID: ${batch_id}\n")
        f.write(f"Cells (raw): {adata.n_obs}\n")
        f.write(f"Genes (features): {adata.n_vars}\n")
        f.write(f"Median UMI per cell: {adata.obs['n_counts'].median():.0f}\n")
        f.write(f"Median genes per cell: {adata.obs['n_genes'].median():.0f}\n")
        f.write(f"Min UMI: {adata.obs['n_counts'].min():.0f}\n")
        f.write(f"Max UMI: {adata.obs['n_counts'].max():.0f}\n")
    
    print(f"Aggregated: {adata.n_obs} cells × {adata.n_vars} genes")
    print(f"Mean UMI per cell: {adata.obs['n_counts'].mean():.0f}")
    PYSCRIPT
    """
}
