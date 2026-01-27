// modules/metadata.nf
// HCA metadata compliance: validation, augmentation, and export

process METADATA_VALIDATE {
    label 'cpu_small'
    publishDir "${params.outdir}/metadata/validation", mode: 'copy'

    input:
    tuple val(batch_id), path(filtered_h5ad)

    output:
    tuple val(batch_id), path(filtered_h5ad), emit: validated
    path("metadata_validation.txt"), emit: report
    path("metadata_metrics.csv"), emit: metrics

    script:
    """
    python3 << 'PYSCRIPT'
    import logging
    import scanpy as sc
    import pandas as pd
    from datetime import datetime
    import warnings
    warnings.filterwarnings('ignore')

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger('METADATA_VALIDATE')

    logger.info("=" * 70)
    logger.info("HCA Metadata Validation")
    logger.info("=" * 70)

    start_time = datetime.now()

    adata = sc.read_h5ad("${filtered_h5ad}")
    logger.info(f"Loaded {adata.n_obs} cells for validation")

    required = ['organism', 'donor_id', 'tissue', 'batch']
    recommended = ['sample_id', 'assay', 'library_preparation_protocol']

    rows = []
    for field in required + recommended:
        if field in adata.obs:
            total = adata.n_obs
            missing = int(adata.obs[field].isna().sum())
            pct_missing = missing / total if total else 0
            n_unique = adata.obs[field].nunique(dropna=True)
            rows.append({
                'field': field,
                'status': 'present',
                'pct_missing': pct_missing,
                'n_unique': int(n_unique)
            })
        else:
            rows.append({
                'field': field,
                'status': 'missing',
                'pct_missing': 1.0,
                'n_unique': 0
            })

    metrics_df = pd.DataFrame(rows)
    metrics_df.to_csv("metadata_metrics.csv", index=False)

    with open("metadata_validation.txt", "w") as f:
        f.write("HCA Metadata Validation Report\n")
        f.write("=" * 70 + "\n")
        f.write(f"Batch: ${batch_id}\n")
        f.write(f"Cells: {adata.n_obs}\n\n")
        f.write("Field Status:\n")
        for _, row in metrics_df.iterrows():
            f.write(f"  {row['field']}: {row['status']} (missing={row['pct_missing']:.2%}, unique={row['n_unique']})\n")

    elapsed = (datetime.now() - start_time).total_seconds()
    logger.info(f"Validation finished in {elapsed:.1f}s")
    logger.info("=" * 70)
    PYSCRIPT
    """
}

process METADATA_INJECT {
    label 'cpu_medium'
    publishDir "${params.outdir}/analysis/final", mode: 'copy'

    input:
    tuple val(batch_id), path(clustered_h5ad)
    val(all_samples)

    output:
    tuple val(batch_id), path("annotated.h5ad"), emit: annotated
    path("metadata_injection_report.txt"), emit: report
    path("metadata_injection_metrics.csv"), emit: metrics

    script:
    """
    python3 << 'PYSCRIPT'
    import logging
    import scanpy as sc
    import pandas as pd
    import json
    from datetime import datetime
    import warnings
    warnings.filterwarnings('ignore')

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger('METADATA_INJECT')

    logger.info("=" * 70)
    logger.info("HCA Metadata Injection")
    logger.info("=" * 70)

    start_time = datetime.now()

    adata = sc.read_h5ad("${clustered_h5ad}")
    logger.info(f"Loaded {adata.n_obs} cells for annotation")

    # all_samples is expected to be a JSON-serializable list of dicts
    try:
        samples_list = json.loads('''${all_samples}''')
    except Exception:
        logger.warning("Failed to parse all_samples as JSON; using empty list")
        samples_list = []

    sample_lookup = {}
    for sample in samples_list:
        if sample.get('batch_id') == '${batch_id}':
            sample_lookup[sample.get('sample_id', '')] = sample

    # Initialize required HCA fields
    adata.obs['organism'] = "${params.organism}"
    adata.obs['tissue'] = "${params.organ}"
    adata.obs['batch'] = "${batch_id}"
    if 'donor_id' not in adata.obs:
        adata.obs['donor_id'] = 'unknown'
    if 'sample_id' not in adata.obs:
        adata.obs['sample_id'] = 'unknown'

    def assign_sample(cell_barcode: str) -> str:
        for sample_id in sample_lookup.keys():
            if cell_barcode.startswith(sample_id):
                return sample_id
        return 'unknown'

    # Map sample_id by barcode prefix
    adata.obs['sample_id'] = adata.obs_names.to_series().apply(assign_sample)

    # Fill donor/tissue from sample metadata when available
    for sample_id, meta in sample_lookup.items():
        mask = adata.obs['sample_id'] == sample_id
        if mask.any():
            if 'donor_id' in meta:
                adata.obs.loc[mask, 'donor_id'] = meta['donor_id']
            if 'tissue' in meta:
                adata.obs.loc[mask, 'tissue'] = meta['tissue']
            if 'assay' in meta:
                adata.obs.loc[mask, 'assay'] = meta['assay']
            if 'library_preparation_protocol' in meta:
                adata.obs.loc[mask, 'library_preparation_protocol'] = meta['library_preparation_protocol']

    # Processing metadata
    adata.uns['processing'] = {
        'timestamp': start_time.isoformat(),
        'pipeline': 'MAPseq-nextflow',
        'version': '1.0.0',
        'batch_id': '${batch_id}',
        'software_versions': {
            'scanpy': sc.__version__,
            'cellranger': '8.0.0'
        }
    }

    adata.uns['hca_compliance'] = {
        'format': 'anndata_h5ad',
        'organism': '${params.organism}',
        'tissue': '${params.organ}',
        'processing_level': 'integrated',
        'batch_id': '${batch_id}'
    }

    # Save annotated object
    adata.write_h5ad("annotated.h5ad", compression="gzip")

    # Metrics
    field_presence = {field: int(field in adata.obs) for field in ['organism', 'donor_id', 'tissue', 'batch', 'sample_id']}
    metrics = {
        'batch_id': '${batch_id}',
        'cells': int(adata.n_obs),
        **field_presence
    }
    pd.DataFrame([metrics]).to_csv("metadata_injection_metrics.csv", index=False)

    # Report
    with open("metadata_injection_report.txt", "w") as f:
        f.write("HCA Metadata Injection Report\n")
        f.write("=" * 70 + "\n")
        f.write(f"Batch: ${batch_id}\n")
        f.write(f"Cells: {adata.n_obs}\n")
        f.write(f"Samples in lookup: {len(sample_lookup)}\n")
        f.write("Fields present: " + ", ".join([k for k, v in field_presence.items() if v]) + "\n")

    elapsed = (datetime.now() - start_time).total_seconds()
    logger.info(f"Metadata injection finished in {elapsed:.1f}s")
    logger.info("=" * 70)
    PYSCRIPT
    """
}
