# Sequencing Portal Integration Guide

This pipeline now emits a portal-friendly manifest per batch for ingestion by the sequencing portal (e.g., `/home/elliott/Apps/github/sequencing_portal`). The manifest captures stable relative paths, checksums, sizes, and labels for all key outputs.

## Where to find the manifest
- Location per batch: `portal/<batch_id>/`
- Files produced:
  - `portal_manifest.json` — primary machine-readable manifest
  - `portal_manifest.csv` — tabular view of the same entries
  - `portal_manifest_report.txt` — summary and missing-file notes

## JSON schema (v1.0)
```json
{
  "schema_version": "1.0",
  "generated_at": "2026-01-27T12:34:56Z",
  "batch_id": "B001",
  "outdir": "/abs/path/to/outdir",
  "entries": [
    {
      "batch_id": "B001",
      "kind": "expression",           // expression | qc | vdj | metadata | metrics | diagnostic | report
      "label": "clustered_h5ad",      // stable label per file type
      "path": "analysis/clustering/B001/clustered.h5ad", // relative to outdir
      "bytes": 123456,
      "md5": "<checksum>",
      "mtime": "2026-01-27T12:34:00Z",
      "source_process": "RAPIDS_CLUSTERING",
      "schema_version": "1.0",
      "tags": ["h5ad", "clustered"],
      "note": "optional"
    }
  ]
}
```

## File coverage (entries are included when files exist)
- Expression / analysis:
  - `analysis/clustering/<batch>/clustered.h5ad`
  - `analysis/harmony/<batch>/harmony_integrated.h5ad`
  - `analysis/final/<batch>/annotated.h5ad`
  - `analysis/final/<batch>/vdj_integrated.h5ad`
- QC:
  - `analysis/qc/<batch>/qc_filtered.h5ad`
  - `analysis/qc/<batch>/doublet_metrics.csv`
  - `analysis/qc/<batch>/doublet_diagnostic.csv`
  - `analysis/qc/<batch>/cellbender_metrics.csv`
- Metrics / reports:
  - `analysis/clustering/<batch>/clustering_metrics.csv`
  - `analysis/clustering/<batch>/clustering_report.txt`
  - `analysis/clustering/<batch>/cluster_summary.csv`
  - `analysis/harmony/<batch>/harmony_metrics.csv`
  - `analysis/harmony/<batch>/harmony_report.txt`
  - `analysis/harmony/<batch>/harmony_diagnostics.csv`
  - `analysis/final/<batch>/metadata_injection_metrics.csv`
  - `analysis/final/<batch>/metadata_injection_report.txt`
- VDJ (when VDJ enabled):
  - `vdj/clonotypes/vdj_clonotypes.csv`
  - `vdj/clonotypes/vdj_metrics.csv`
  - `vdj/clonotypes/vdj_diagnostics.csv`
  - `analysis/final/<batch>/vdj_integration_metrics.csv`
  - `analysis/final/<batch>/vdj_cluster_summary.csv`
  - `analysis/final/<batch>/vdj_integration_report.txt`

## Ingestion tips for the portal
1. Read `portal_manifest.json` and resolve each `path` relative to `outdir`.
2. Use `kind` + `label` to route to the correct component (e.g., H5AD loader vs CSV parser).
3. Verify integrity with `md5` and `bytes` before ingesting.
4. Use `tags` for UI grouping (e.g., `vdj`, `metrics`, `report`).
5. Honor missing entries: only present files are listed. Empty lists mean optional steps were skipped.

## Example Python ingestion snippet
```python
import json
from pathlib import Path

manifest = json.loads(Path('portal/B001/portal_manifest.json').read_text())
outdir = Path(manifest['outdir'])
for entry in manifest['entries']:
    full_path = outdir / entry['path']
    print(entry['label'], full_path)
```

## How to extend
- To add new file types, edit `modules/portal_manifest.nf` and append an `add_entry(...)` call with a stable `label`.
- Keep `schema_version` bumped when changing structure.
- Prefer relative paths inside `outdir` to keep portability across environments.
