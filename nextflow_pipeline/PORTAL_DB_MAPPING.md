# Portal Database Mapping (Guidance)

Path: database schema provided at ../MAPseq/database_schema.sql (outside pipeline repo); key tables: investigations, studies, assays, samples, assay_samples, metrics, files.

## Minimal ingestion flow
1) Resolve portal_manifest.json (per batch) → get outdir + entries.
2) Create or reuse hierarchy (investigation → study → assay) according to your portal rules; example defaults:
   - investigation.name = project or run label
   - study.name = genome or experiment label
   - assay.name = "batch_<batch_id>"
   - assays.assay_type = 'scRNA-seq', library_type = 'GEX', vdj_type = 'BCR/TCR' when VDJ present
3) For each sample/barcode group:
   - samples.name = sample_id (or sample prefix) from pipeline metadata
   - assay_samples rows link assay ↔ sample; store barcode prefix and set quality_status ('pending' by default)
4) Files table: insert one row per manifest entry
   - file_path = manifest entry path (relative to outdir or absolute)
   - file_type: map from `kind`/`label` (see table below)
   - file_size = bytes; checksum = md5
   - metrics_file_path: you may store main metrics CSV/JSON reference in assay_samples if desired
5) Metrics table: parse select CSVs and insert numeric/enum metrics with metric_category tags.

## File-type mapping (manifest → files.file_type)
- clustered_h5ad → expression_clustered
- harmony_integrated_h5ad → expression_integrated
- annotated_h5ad → expression_annotated
- vdj_integrated_h5ad → vdj_annotated
- qc_filtered_h5ad → qc_filtered
- doublet_metrics.csv → qc_metrics_doublet
- doublet_diagnostic.csv → qc_diag_doublet
- cellbender_metrics.csv → qc_metrics_cellbender
- clustering_metrics.csv → metrics_clustering
- clustering_report.txt → report_clustering
- cluster_summary.csv → summary_cluster
- harmony_metrics.csv → metrics_harmony
- harmony_report.txt → report_harmony
- harmony_diagnostics.csv → diag_harmony
- metadata_injection_metrics.csv → metrics_metadata
- metadata_injection_report.txt → report_metadata
- vdj_clonotypes.csv → vdj_clonotypes
- vdj_metrics.csv → metrics_vdj
- vdj_diagnostics.csv → diag_vdj
- vdj_integration_metrics.csv → metrics_vdj_integration
- vdj_cluster_summary.csv → summary_vdj_cluster
- vdj_integration_report.txt → report_vdj

## Metrics suggestions
Insert into metrics table (assay_sample scoped where possible; otherwise assay-level with a synthetic sample like "batch-<id>"):
- From clustering_metrics.csv: n_cells, n_genes, n_clusters, leiden_resolution
- From harmony_metrics.csv: n_cells, n_batches, n_hvg, harmony_status, harmony_pcs
- From vdj_integration_metrics.csv: cells_with_vdj, productive_cells, clonotypes
- From doublet_metrics.csv: doublet_rate, cells_kept
- From cellbender_metrics.csv: ambient_percent_before/after, cells_retained

## Join keys / identifiers
- batch_id: use as assay.name or suffix
- sample_id: derived from barcode prefixes (pipeline metadata injects sample_id when possible)
- barcode: `vdj_barcode` / obs_names stripped of `-1`; use for assay_samples.barcode when populating per-barcode rows (optional)

## Ingestion pseudocode (Python snippet)
```python
from pathlib import Path
import json, csv

manifest = json.loads(Path('portal/B001/portal_manifest.json').read_text())
outdir = Path(manifest['outdir'])
entries = manifest['entries']

files_rows = []
for e in entries:
    files_rows.append({
        'file_path': str(outdir / e['path']),
        'file_type': e['label'],
        'file_size': e['bytes'],
        'checksum': e['md5'],
    })
# Insert files_rows into portal DB (files table)
```

## Actions needed in portal code
- Add a manifest loader that reads portal_manifest.json and resolves paths relative to outdir.
- Map `label` → `file_type` using the table above.
- Parse key metrics CSVs and upsert into metrics table.
- Optionally create an assay per batch and link samples via assay_samples.

## Notes
- Manifest uses relative paths under `params.outdir`; keep that stable for portability.
- If schema evolves, bump manifest schema_version and adjust mapping.
