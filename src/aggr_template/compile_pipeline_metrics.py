#!/usr/bin/env python3
"""Compile MAPseq mapping, demultiplexing, and QC CSVs into a long-form ledger."""

import argparse
import csv
from pathlib import Path


OUTPUT_FIELDS = [
    "project_id",
    "run_id",
    "modality",
    "stage",
    "sample_id",
    "sample_scope",
    "metric",
    "value",
    "unit",
    "source_file",
]
INPUT_PATTERNS = (
    "mapping.metrics_*.csv",
    "demux_metrics.*.csv",
    "library_stats.*.csv",
    "QC.*_sampleID_filtering.stats.csv",
)
IDENTIFIER_COLUMNS = {
    "run_id",
    "sample_id",
    "library_id",
    "atac_id",
    "asap_id",
    "project_id",
    "modality",
}


def normalize_header(value):
    return "_".join(value.strip().lower().replace(" ", "_").split())


def infer_stage(path):
    name = path.name
    if name.startswith("mapping.metrics_"):
        return "mapping"
    if name.startswith("QC."):
        return "qc"
    return "demultiplexing"


def infer_modality(path):
    text = "/".join(path.parts).lower()
    if "atac" in text or "asap" in text:
        return "ATAC"
    if "rna" in text or "vdj" in text:
        return "RNA"
    return "unknown"


def read_records(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return [], []
        headers = [normalize_header(name or "") for name in reader.fieldnames]
        records = []
        for raw_record in reader:
            records.append(
                {
                    normalize_header(key or ""): value
                    for key, value in raw_record.items()
                    if key is not None
                }
            )
    return headers, records


def choose_sample_id(record, fallback):
    for column in ("sample_id", "library_id", "atac_id"):
        value = record.get(column, "")
        if value:
            return value
    return fallback


def rows_from_file(path, root, project_id):
    headers, records = read_records(path)
    if not headers:
        return []

    relative_path = path.relative_to(root)
    fallback_sample = ""
    if path.name.startswith("mapping.metrics_"):
        fallback_sample = path.stem.removeprefix("mapping.metrics_")
    inferred_run = relative_path.parts[0] if len(relative_path.parts) > 1 else ""
    output = []
    metric_name_column = next(
        (name for name in ("metric", "metric_name") if name in headers), None
    )
    metric_value_column = next(
        (name for name in ("value", "metric_value") if name in headers), None
    )

    for record in records:
        sample_id = choose_sample_id(record, fallback_sample)
        if not sample_id:
            continue
        run_id = record.get("run_id") or inferred_run
        common = {
            "project_id": project_id,
            "run_id": run_id,
            "modality": record.get("modality") or infer_modality(relative_path),
            "stage": infer_stage(path),
            "sample_id": sample_id,
            "sample_scope": {
                "mapping": "cellranger_sample",
                "demultiplexing": "library",
                "qc": "qc_sample",
            }[infer_stage(path)],
            "unit": "",
            "source_file": relative_path.as_posix(),
        }

        if metric_name_column and metric_value_column:
            metric = record.get(metric_name_column, "")
            value = record.get(metric_value_column, "")
            if metric and value != "":
                output.append({**common, "metric": metric, "value": value})
            continue

        for column in headers:
            if column in IDENTIFIER_COLUMNS or not column:
                continue
            value = record.get(column, "")
            if value == "":
                continue
            output.append({**common, "metric": column, "value": value})
    return output


def compile_metrics(root, output_path, project_id):
    root = root.resolve()
    output_path = output_path.resolve()
    source_files = sorted(
        {
            path
            for pattern in INPUT_PATTERNS
            for path in root.rglob(pattern)
            if path.resolve() != output_path
        }
    )
    rows = []
    for path in source_files:
        rows.extend(rows_from_file(path, root, project_id))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    return len(source_files), len(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "root",
        type=Path,
        help="Common directory containing the run, aggregation, and analysis outputs",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output CSV (default: ROOT/reports/sample_metrics.csv)",
    )
    parser.add_argument("--project-id", help="Override the project ID in output rows")
    args = parser.parse_args()

    root = args.root.resolve()
    if not root.is_dir():
        parser.error(f"not a directory: {root}")
    output_path = args.output or root / "reports" / "sample_metrics.csv"
    file_count, row_count = compile_metrics(
        root, output_path, args.project_id or root.name
    )
    print(f"Wrote {row_count} metric rows from {file_count} source CSVs to {output_path}")


if __name__ == "__main__":
    main()