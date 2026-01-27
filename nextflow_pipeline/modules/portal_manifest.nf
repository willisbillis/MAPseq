// modules/portal_manifest.nf
// Generate portal-ready manifest (JSON + CSV) for downstream ingestion

process PORTAL_MANIFEST {
    label 'cpu_small'
    publishDir "${params.outdir}/portal/${batch_id}", mode: 'copy'

    input:
    tuple val(batch_id), path(final_h5ad)

    output:
    path("portal_manifest.json"), emit: manifest_json
    path("portal_manifest.csv"), emit: manifest_csv
    path("portal_manifest_report.txt"), emit: manifest_report

    script:
    """
    python3 << 'PYSCRIPT'
    import json
    import hashlib
    from pathlib import Path
    from datetime import datetime
    import pandas as pd
    import logging
    import os

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger('PORTAL_MANIFEST')

    portal_schema_version = "1.0"
    outdir = Path("${params.outdir}").resolve()
    batch_id = "${batch_id}"
    run_ts = datetime.now().isoformat()

    def md5sum(path: Path, chunk=8192):
        h = hashlib.md5()
        with path.open('rb') as f:
            for block in iter(lambda: f.read(chunk), b''):
                h.update(block)
        return h.hexdigest()

    def add_entry(entries, path: Path, kind: str, label: str, source: str, tags=None, extra=None):
        if not path.exists():
            return
        tags = tags or []
        extra = extra or {}
        rel = str(path.relative_to(outdir)) if path.is_absolute() and outdir in path.parents else str(path)
        stat = path.stat()
        entries.append({
            'batch_id': batch_id,
            'kind': kind,
            'label': label,
            'path': rel,
            'bytes': stat.st_size,
            'md5': md5sum(path),
            'mtime': datetime.fromtimestamp(stat.st_mtime).isoformat(),
            'source_process': source,
            'schema_version': portal_schema_version,
            'tags': tags,
            **extra
        })

    entries = []

    # Core analysis outputs
    clustering_dir = outdir / "analysis" / "clustering" / batch_id
    harmony_dir = outdir / "analysis" / "harmony" / batch_id
    final_dir = outdir / "analysis" / "final" / batch_id
    qc_dir = outdir / "analysis" / "qc" / batch_id
    vdj_dir = outdir / "vdj" / "clonotypes"

    add_entry(entries, clustering_dir / "clustered.h5ad", "expression", "clustered_h5ad", "RAPIDS_CLUSTERING", tags=["h5ad", "clustered"])
    add_entry(entries, harmony_dir / "harmony_integrated.h5ad", "expression", "harmony_integrated_h5ad", "HARMONY_INTEGRATION", tags=["h5ad", "integrated"])
    add_entry(entries, final_dir / "annotated.h5ad", "metadata", "annotated_h5ad", "METADATA_INJECT", tags=["h5ad", "annotated"])
    add_entry(entries, final_dir / "vdj_integrated.h5ad", "vdj", "vdj_integrated_h5ad", "VDJ_INTEGRATION", tags=["h5ad", "vdj"])

    # QC artifacts
    add_entry(entries, qc_dir / "qc_filtered.h5ad", "qc", "qc_filtered_h5ad", "SCDBLFINDER_FILTER", tags=["h5ad", "qc"])
    add_entry(entries, qc_dir / "doublet_metrics.csv", "qc", "doublet_metrics", "SCDBLFINDER_FILTER", tags=["metrics", "csv"])
    add_entry(entries, qc_dir / "doublet_diagnostic.csv", "qc", "doublet_diagnostic", "SCDBLFINDER_FILTER", tags=["diagnostic", "csv"])
    add_entry(entries, qc_dir / "cellbender_metrics.csv", "qc", "cellbender_metrics", "CELLBENDER_QC", tags=["metrics", "csv"])

    # Analysis metrics
    add_entry(entries, clustering_dir / "clustering_metrics.csv", "metrics", "clustering_metrics", "RAPIDS_CLUSTERING", tags=["metrics", "csv"])
    add_entry(entries, clustering_dir / "clustering_report.txt", "report", "clustering_report", "RAPIDS_CLUSTERING", tags=["report"])
    add_entry(entries, clustering_dir / "cluster_summary.csv", "metrics", "cluster_summary", "RAPIDS_CLUSTERING", tags=["summary", "csv"])

    add_entry(entries, harmony_dir / "harmony_metrics.csv", "metrics", "harmony_metrics", "HARMONY_INTEGRATION", tags=["metrics", "csv"])
    add_entry(entries, harmony_dir / "harmony_report.txt", "report", "harmony_report", "HARMONY_INTEGRATION", tags=["report"])
    add_entry(entries, harmony_dir / "harmony_diagnostics.csv", "diagnostic", "harmony_diagnostics", "HARMONY_INTEGRATION", tags=["diagnostic", "csv"])

    # Metadata artifacts
    add_entry(entries, final_dir / "metadata_injection_metrics.csv", "metadata", "metadata_metrics", "METADATA_INJECT", tags=["metrics", "csv"])
    add_entry(entries, final_dir / "metadata_injection_report.txt", "report", "metadata_report", "METADATA_INJECT", tags=["report"])

    # VDJ artifacts
    add_entry(entries, vdj_dir / "vdj_clonotypes.csv", "vdj", "vdj_clonotypes", "VDJ_CLONOTYPE_ANALYSIS", tags=["vdj", "csv"])
    add_entry(entries, vdj_dir / "vdj_metrics.csv", "vdj", "vdj_metrics", "VDJ_CLONOTYPE_ANALYSIS", tags=["metrics", "csv"])
    add_entry(entries, vdj_dir / "vdj_diagnostics.csv", "vdj", "vdj_diagnostics", "VDJ_CLONOTYPE_ANALYSIS", tags=["diagnostic", "csv"])
    add_entry(entries, final_dir / "vdj_integration_metrics.csv", "vdj", "vdj_integration_metrics", "VDJ_INTEGRATION", tags=["metrics", "csv"])
    add_entry(entries, final_dir / "vdj_cluster_summary.csv", "vdj", "vdj_cluster_summary", "VDJ_INTEGRATION", tags=["summary", "csv"])
    add_entry(entries, final_dir / "vdj_integration_report.txt", "report", "vdj_integration_report", "VDJ_INTEGRATION", tags=["report"])

    # Final input pointer (as received by this process)
    staged_final = Path("${final_h5ad}").resolve()
    add_entry(entries, staged_final, "expression", "staged_final_input", "PORTAL_MANIFEST", tags=["h5ad", "staged"], extra={'note': 'staged workdir copy'})

    manifest = {
        'schema_version': portal_schema_version,
        'generated_at': run_ts,
        'batch_id': batch_id,
        'outdir': str(outdir),
        'entries': entries
    }

    Path("portal_manifest.json").write_text(json.dumps(manifest, indent=2))
    pd.DataFrame(entries).to_csv("portal_manifest.csv", index=False)

    missing = [e['label'] for e in entries if not (outdir / e['path']).exists()]
    present = len(entries) - len(missing)

    with open("portal_manifest_report.txt", "w") as f:
        f.write("Portal Manifest Report\n")
        f.write("=" * 70 + "\n")
        f.write(f"Batch ID: {batch_id}\n")
        f.write(f"Outdir: {outdir}\n")
        f.write(f"Entries written: {len(entries)}\n")
        f.write(f"Missing files (skipped): {len(missing)}\n")
        if missing:
            f.write("Missing labels: " + ", ".join(missing) + "\n")

    logger.info(f"Manifest entries: {len(entries)} (missing counted in report)")
    logger.info("Portal manifest generation complete")
    PYSCRIPT
    """
}
