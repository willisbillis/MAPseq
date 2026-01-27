// modules/vdj.nf
// Production-grade V(D)J processing with scirpy (if available) and robust fallbacks

process VDJ_CLONOTYPE_ANALYSIS {
    label 'cpu_medium'
    publishDir "${params.outdir}/vdj/clonotypes", mode: 'copy'

    input:
    path(vdj_outputs)

    output:
    path("vdj_clonotypes.csv"), emit: clonotypes
    path("vdj_report.txt"), emit: report
    path("vdj_metrics.csv"), emit: metrics
    path("vdj_diagnostics.csv"), emit: diagnostics

    script:
    """
    python3 << 'PYSCRIPT'
    import logging
    from pathlib import Path
    import pandas as pd
    import numpy as np
    from datetime import datetime
    import warnings
    warnings.filterwarnings('ignore')

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger('VDJ_CLONOTYPE_ANALYSIS')

    logger.info("=" * 70)
    logger.info("V(D)J Clonotype Analysis")
    logger.info("=" * 70)

    start_time = datetime.now()

    vdj_dir = Path("${vdj_outputs}")
    contig_files = sorted(vdj_dir.glob("**/filtered_contig*.csv"))
    logger.info(f"Found {len(contig_files)} VDJ contig files")

    metrics = {
        'n_contig_files': len(contig_files),
        'n_contigs': 0,
        'n_cells': 0,
        'n_productive': 0,
        'n_clonotypes': 0,
        'method': 'scirpy'  # updated later if fallback
    }

    if not contig_files:
        logger.warning("No VDJ data found; emitting empty outputs")
        Path("vdj_clonotypes.csv").write_text("barcode,clonotype_id,chain,is_productive,umi_count,cdr3s_aa\n")
        Path("vdj_diagnostics.csv").write_text("metric,value\nno_vdj_data,true\n")
        pd.DataFrame([metrics]).to_csv("vdj_metrics.csv", index=False)
        Path("vdj_report.txt").write_text("No VDJ data to process\n")
    else:
        # Try scirpy if available, otherwise fallback to pandas-based assembly
        use_scirpy = False
        clonotype_df = None
        diag_rows = []

        try:
            import scirpy as ir
            logger.info("Using scirpy for VDJ processing")
            adata_ir = ir.io.read_10x_vdj([str(f) for f in contig_files])
            logger.info(f"Loaded {adata_ir.n_obs} contigs across {adata_ir.obs['barcode'].nunique()} cells")

            # Chain QC and basic filtering
            ir.pp.chain_qc(adata_ir)
            ir.tl.chain_qc(adata_ir)

            # Assign clonotypes
            ir.tl.clonotypes(adata_ir)
            clonotype_df = ir.get.tcr_clonotypes(adata_ir, key='clonotype')

            # scirpy returns index as clonotype_id; reset for CSV
            clonotype_df = clonotype_df.reset_index().rename(columns={'index': 'clonotype_id'})
            clonotype_df.rename(columns={'cell_id': 'barcode'}, inplace=True)

            use_scirpy = True
            metrics['method'] = 'scirpy'
            metrics['n_contigs'] = adata_ir.n_obs
            metrics['n_cells'] = adata_ir.obs['barcode'].nunique()
            metrics['n_productive'] = int((adata_ir.obs['productive'] == 'True').sum()) if 'productive' in adata_ir.obs else 0
            metrics['n_clonotypes'] = clonotype_df['clonotype_id'].nunique()

            # Diagnostics: chain usage
            if 'chain_pairing' in adata_ir.obs:
                chain_counts = adata_ir.obs['chain_pairing'].value_counts()
                for chain, count in chain_counts.items():
                    diag_rows.append({'metric': f"chain_pairing_{chain}", 'value': count})

        except ImportError:
            logger.warning("scirpy not installed; using pandas fallback")
        except Exception as e:
            logger.error(f"scirpy processing failed: {str(e)}", exc_info=True)

        if clonotype_df is None:
            metrics['method'] = 'pandas_fallback'
            logger.info("Running pandas-based fallback assembly")
            all_contigs = []
            for contig_file in contig_files:
                try:
                    df = pd.read_csv(contig_file)
                    all_contigs.append(df)
                except Exception as e:
                    logger.warning(f"Failed to read {contig_file}: {e}")
            if all_contigs:
                contigs_df = pd.concat(all_contigs, ignore_index=True)
                metrics['n_contigs'] = len(contigs_df)
                contigs_df['barcode_std'] = contigs_df['barcode'].str.replace('-1', '', regex=False)

                agg_cols = {
                    'clonotype_id': 'first',
                    'chain': lambda x: ';'.join(sorted(set(x.dropna()))),
                    'productive': lambda x: any(x.astype(str).str.lower() == 'true'),
                    'umi_count': lambda x: x.fillna(0).astype(int).sum() if 'umi_count' in contigs_df else x.count(),
                    'cdr3': lambda x: ';'.join(sorted(set(x.dropna()))),
                    'cdr3_nt': lambda x: ';'.join(sorted(set(x.dropna())))
                }

                clonotype_df = contigs_df.groupby('barcode_std').agg(agg_cols).reset_index().rename(columns={'barcode_std': 'barcode'})
                clonotype_df['is_productive'] = clonotype_df['productive'].astype(bool)
                clonotype_df['umi_count'] = clonotype_df.get('umi_count', 0).astype(int)
                clonotype_df['cdr3s_aa'] = clonotype_df.get('cdr3', '')

                metrics['n_cells'] = clonotype_df['barcode'].nunique()
                metrics['n_productive'] = int(clonotype_df['is_productive'].sum())
                metrics['n_clonotypes'] = clonotype_df['clonotype_id'].nunique()

                diag_rows.append({'metric': 'fallback_used', 'value': True})
            else:
                logger.warning("No contigs could be read; emitting empty outputs")
                clonotype_df = pd.DataFrame(columns=['barcode', 'clonotype_id', 'chain', 'is_productive', 'umi_count', 'cdr3s_aa'])

        # Finalize outputs
        clonotype_df = clonotype_df[['barcode', 'clonotype_id', 'chain', 'is_productive', 'umi_count', 'cdr3s_aa']].copy()
        clonotype_df['clonotype_id'] = clonotype_df['clonotype_id'].fillna('unassigned')
        clonotype_df['chain'] = clonotype_df['chain'].fillna('unknown')
        clonotype_df['is_productive'] = clonotype_df['is_productive'].fillna(False).astype(bool)
        clonotype_df['umi_count'] = clonotype_df['umi_count'].fillna(0).astype(int)
        clonotype_df['cdr3s_aa'] = clonotype_df['cdr3s_aa'].fillna('')

        clonotype_df.to_csv("vdj_clonotypes.csv", index=False)

        # Diagnostics
        diag_df = pd.DataFrame(diag_rows) if diag_rows else pd.DataFrame({'metric': ['analysis_method'], 'value': [metrics['method']]})
        diag_df.to_csv("vdj_diagnostics.csv", index=False)

        # Metrics
        pd.DataFrame([metrics]).to_csv("vdj_metrics.csv", index=False)

        # Report
        with open("vdj_report.txt", "w") as f:
            f.write("VDJ Clonotype Analysis Report\n")
            f.write("=" * 70 + "\n")
            f.write(f"Input directory: {vdj_dir}\n")
            f.write(f"Contig files: {metrics['n_contig_files']}\n")
            f.write(f"Total contigs: {metrics['n_contigs']}\n")
            f.write(f"Cells with contigs: {metrics['n_cells']}\n")
            f.write(f"Productive cells: {metrics['n_productive']}\n")
            f.write(f"Unique clonotypes: {metrics['n_clonotypes']}\n")
            f.write(f"Method: {metrics['method']}\n")

    elapsed = (datetime.now() - start_time).total_seconds()
    logger.info(f"VDJ processing finished in {elapsed:.1f}s")
    logger.info("=" * 70)
    PYSCRIPT
    """
}

process VDJ_INTEGRATION {
    label 'cpu_medium'
    publishDir "${params.outdir}/analysis/final", mode: 'copy'

    input:
    tuple val(batch_id), path(clustered_h5ad)
    path(vdj_clonotypes)

    output:
    tuple val(batch_id), path("vdj_integrated.h5ad"), emit: integrated
    path("vdj_integration_report.txt"), emit: report
    path("vdj_integration_metrics.csv"), emit: metrics
    path("vdj_cluster_summary.csv"), emit: summary

    script:
    """
    python3 << 'PYSCRIPT'
    import logging
    import pandas as pd
    import scanpy as sc
    from datetime import datetime
    import numpy as np
    import warnings
    warnings.filterwarnings('ignore')

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger('VDJ_INTEGRATION')

    logger.info("=" * 70)
    logger.info("V(D)J Integration into primary analysis")
    logger.info("=" * 70)

    start_time = datetime.now()

    # Load data
    adata = sc.read_h5ad("${clustered_h5ad}")
    clonotypes = pd.read_csv("${vdj_clonotypes}")

    logger.info(f"Cells in expression object: {adata.n_obs}")
    logger.info(f"Clonotype rows: {len(clonotypes)}")

    def standardize_barcode(bc: str) -> str:
        if pd.isna(bc):
            return ''
        bc = str(bc)
        if bc.endswith('-1'):
            bc = bc[:-2]
        return bc

    adata.obs['vdj_barcode'] = [standardize_barcode(bc) for bc in adata.obs_names]
    clonotypes['barcode_std'] = clonotypes['barcode'].apply(standardize_barcode)

    merged_obs = adata.obs.merge(
        clonotypes,
        left_on='vdj_barcode',
        right_on='barcode_std',
        how='left'
    )

    merged_obs['clonotype_id'] = merged_obs['clonotype_id'].fillna('uncovered')
    merged_obs['is_productive'] = merged_obs['is_productive'].fillna(False).astype(bool)
    merged_obs['umi_count'] = merged_obs['umi_count'].fillna(0).astype(int)
    merged_obs['cdr3s_aa'] = merged_obs['cdr3s_aa'].fillna('')
    merged_obs['vdj_present'] = merged_obs['clonotype_id'] != 'uncovered'

    adata.obs = merged_obs

    n_cells_with_vdj = int(adata.obs['vdj_present'].sum())
    n_clonotypes = int(adata.obs['clonotype_id'].nunique())
    n_productive = int(adata.obs['is_productive'].sum())

    adata.uns['vdj_summary'] = {
        'n_cells': int(adata.n_obs),
        'n_cells_with_vdj': n_cells_with_vdj,
        'n_productive_cells': n_productive,
        'n_clonotypes': n_clonotypes,
        'timestamp': start_time.isoformat(),
        'barcode_mapping': 'stripped -1 suffix for matching'
    }

    # Cluster summary if clustering available
    cluster_key = 'leiden' if 'leiden' in adata.obs else None
    cluster_rows = []
    if cluster_key:
        logger.info(f"Summarizing VDJ per cluster using {cluster_key}")
        cluster_counts = adata.obs.groupby(cluster_key).agg({
            'vdj_present': 'sum',
            'is_productive': 'sum',
            'clonotype_id': 'nunique',
            cluster_key: 'size'
        }).rename(columns={cluster_key: 'n_cells'}).reset_index()
        cluster_counts['vdj_fraction'] = cluster_counts['vdj_present'] / cluster_counts['n_cells']
        cluster_counts.to_csv("vdj_cluster_summary.csv", index=False)
        cluster_rows = cluster_counts.to_dict(orient='records')
    else:
        logger.warning("No clustering key found; writing empty cluster summary")
        pd.DataFrame(columns=['cluster', 'n_cells', 'vdj_present', 'vdj_fraction']).to_csv("vdj_cluster_summary.csv", index=False)

    # Save integrated object
    adata.write_h5ad("vdj_integrated.h5ad", compression="gzip")

    # Metrics
    metrics = {
        'batch_id': '${batch_id}',
        'cells_total': int(adata.n_obs),
        'cells_with_vdj': n_cells_with_vdj,
        'productive_cells': n_productive,
        'clonotypes': n_clonotypes,
        'cluster_key': cluster_key or 'none'
    }
    pd.DataFrame([metrics]).to_csv("vdj_integration_metrics.csv", index=False)

    # Report
    with open("vdj_integration_report.txt", "w") as f:
        f.write("V(D)J Integration Report\n")
        f.write("=" * 70 + "\n")
        f.write(f"Batch ID: ${batch_id}\n")
        f.write(f"Cells: {adata.n_obs}\n")
        f.write(f"Cells with VDJ: {n_cells_with_vdj}\n")
        f.write(f"Productive cells: {n_productive}\n")
        f.write(f"Unique clonotypes: {n_clonotypes}\n")
        f.write(f"Clustering key: {cluster_key if cluster_key else 'none'}\n")

    elapsed = (datetime.now() - start_time).total_seconds()
    logger.info(f"VDJ integration finished in {elapsed:.1f}s")
    logger.info("=" * 70)
    PYSCRIPT
    """
}
