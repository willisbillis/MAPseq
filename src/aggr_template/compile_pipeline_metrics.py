#!/usr/bin/env python3
"""Compile MAPseq mapping, demultiplexing, and QC CSVs into a long-form ledger."""

import argparse
import csv
import os
import re
import tempfile
import time
from contextlib import contextmanager
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
LOCK_TIMEOUT_SECONDS = 300
LOCK_POLL_SECONDS = 0.1
QC_METRIC_UNITS = {
    "RNA": {
        "unfiltered_cells": "cells",
        "filtered_cells": "cells",
        "doublet_cells": "cells",
        "negative_cells": "cells",
        "mean_depth_before": "counts_per_cell",
        "mean_depth_after": "counts_per_cell",
        "mean_features_after": "features_per_cell",
        "median_percent_mt_after": "percent",
    },
    "ATAC": {
        "unfiltered_cells": "cells",
        "filtered_cells": "cells",
        "doublet_cells": "cells",
        "negative_cells": "cells",
        "mean_depth_before": "fragments_per_cell",
        "mean_depth_after": "fragments_per_cell",
        "mean_tss_after": "tss_enrichment",
        "mean_pct_reads_in_peaks_after": "percent",
        "mean_nucleosome_signal_after": "score",
    },
}
NUMERIC_VALUE = re.compile(r"-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?%?\Z")


@contextmanager
def metrics_lock(output_path):
    """Serialize ledger writers using the same atomic directory lock as R."""
    lock_path = Path(f"{output_path}.lock")
    deadline = time.monotonic() + LOCK_TIMEOUT_SECONDS
    while True:
        try:
            os.mkdir(lock_path)
            break
        except FileExistsError:
            if time.monotonic() >= deadline:
                raise TimeoutError(
                    f"Timed out waiting for metrics lock: {lock_path}. "
                    "Confirm no writer is running, then remove a stale lock directory if needed."
                )
            time.sleep(LOCK_POLL_SECONDS)

    try:
        yield
    finally:
        os.rmdir(lock_path)


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
    if path.name.lower().startswith("qc.atac") or any(
        part in {"atac", "atac.asap", "atac_demuxing"} for part in path.parts
    ):
        return "ATAC"
    if path.name.lower().startswith("qc.rna") or any(
        part in {"rna", "rna.fb.vdj", "rna_demuxing"} for part in path.parts
    ):
        return "RNA"
    return "unknown"


def safe_source(path, root):
    """Reject symlink inputs and paths resolving outside the scan root."""
    path = Path(path)
    root = Path(root).resolve()
    current = path
    while current != root and current != current.parent:
        if current.is_symlink():
            return False
        current = current.parent
    try:
        path.resolve().relative_to(root)
    except ValueError:
        return False
    return path.exists()


def read_aggregation_manifest(manifest_path, root, modality):
    entries = []
    if not safe_source(manifest_path, root):
        return entries
    with manifest_path.open(newline="", encoding="utf-8-sig") as handle:
        for record in csv.DictReader(handle):
            source_column = "molecule_h5" if modality == "RNA" else "fragments"
            sample_column = "sample_id" if modality == "RNA" else "library_id"
            source_value = (record.get(source_column) or "").strip()
            sample_id = (record.get(sample_column) or "").strip()
            if not source_value or not sample_id:
                continue
            source_path = Path(source_value)
            if not source_path.is_absolute():
                source_path = root / source_path
            if not safe_source(source_path, root):
                continue
            pipeline_dir = next(
                (parent for parent in source_path.parents if parent.name == "pipeline"),
                None,
            )
            if pipeline_dir is None:
                continue
            run_dir = pipeline_dir.parent
            if safe_source(run_dir, root):
                entries.append((run_dir.resolve(), modality, sample_id))
    return entries


def discover_sources(root, project_dir):
    sources = []
    manifests = (
        (project_dir / "pipeline/RNA.FB.VDJ/aggr.csv", "RNA"),
        (project_dir / "pipeline/ATAC.ASAP/ATAC/aggr.csv", "ATAC"),
    )
    for manifest_path, modality in manifests:
        for run_dir, source_modality, sample_id in read_aggregation_manifest(
            manifest_path, root, modality
        ):
            if source_modality == "RNA":
                report_dirs = (run_dir / "pipeline/RNA.FB.VDJ/reports",)
            else:
                report_dirs = (
                    run_dir / "pipeline/ATAC.ASAP/ATAC/reports",
                    run_dir / "pipeline/ATAC.ASAP/reports",
                )
            filename = f"mapping.metrics_{sample_id}.csv"
            candidate = next(
                (
                    report_dir / filename
                    for report_dir in report_dirs
                    if safe_source(report_dir / filename, root)
                    and (report_dir / filename).is_file()
                ),
                None,
            )
            if candidate is not None:
                sources.append((candidate, run_dir.name))

    for pattern in INPUT_PATTERNS:
        for path in project_dir.rglob(pattern):
            if safe_source(path, root) and path.is_file():
                sources.append((path, ""))
    return sources


def spreadsheet_safe(value, field):
    if value is None:
        return ""
    text = str(value)
    stripped = text.lstrip(" \t\r\n")
    if field == "value" and NUMERIC_VALUE.fullmatch(stripped):
        return text
    if stripped.startswith(("=", "+", "-", "@")):
        return "'" + text
    return text


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


def rows_from_file(path, root, project_id, inferred_run=""):
    headers, records = read_records(path)
    if not headers:
        return []

    relative_path = path.relative_to(root)
    fallback_sample = ""
    if path.name.startswith("mapping.metrics_"):
        fallback_sample = path.stem.removeprefix("mapping.metrics_")
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
        modality = record.get("modality") or infer_modality(relative_path)
        common = {
            "project_id": project_id,
            "run_id": run_id,
            "modality": modality,
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
                common["unit"] = QC_METRIC_UNITS.get(modality.upper(), {}).get(metric, "")
                output.append({**common, "metric": metric, "value": value})
            continue

        for column in headers:
            if column in IDENTIFIER_COLUMNS or not column:
                continue
            value = record.get(column, "")
            if value == "":
                continue
            common["unit"] = QC_METRIC_UNITS.get(modality.upper(), {}).get(column, "")
            output.append({**common, "metric": column, "value": value})
    return output


def compile_metrics(root, output_path, project_id, project_dir):
    root = root.resolve()
    output_path = output_path.resolve()
    project_dir = project_dir.resolve()
    try:
        project_dir.relative_to(root)
    except ValueError as error:
        raise ValueError(f"project directory must be inside scan root: {project_dir}") from error
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with metrics_lock(output_path):
        discovered = discover_sources(root, project_dir)
        source_files = sorted(
            {
                (path.resolve(), inferred_run)
                for path, inferred_run in discovered
                if path.resolve() != output_path and safe_source(path, root)
            },
            key=lambda item: (item[0].as_posix(), item[1]),
        )
        rows = []
        for path, inferred_run in source_files:
            rows.extend(rows_from_file(path, root, project_id, inferred_run))

        unique_rows = []
        seen_mapping_rows = set()
        for row in rows:
            if row["stage"] == "mapping":
                identity = tuple(row[field] for field in OUTPUT_FIELDS if field != "source_file")
                if identity in seen_mapping_rows:
                    continue
                seen_mapping_rows.add(identity)
            unique_rows.append(row)
        rows = [
            {field: spreadsheet_safe(row.get(field, ""), field) for field in OUTPUT_FIELDS}
            for row in unique_rows
        ]

        temp_path = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="w",
                newline="",
                encoding="utf-8",
                dir=output_path.parent,
                prefix=f".{output_path.name}.",
                suffix=".tmp",
                delete=False,
            ) as handle:
                temp_path = Path(handle.name)
                writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS)
                writer.writeheader()
                writer.writerows(rows)
            os.replace(temp_path, output_path)
        finally:
            if temp_path is not None:
                temp_path.unlink(missing_ok=True)
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
    parser.add_argument(
        "--project-dir",
        type=Path,
        help="Aggregation project directory; defaults to ROOT/PROJECT_ID",
    )
    parser.add_argument("--project-id", help="Override the project ID in output rows")
    args = parser.parse_args()

    root = args.root.resolve()
    if not root.is_dir():
        parser.error(f"not a directory: {root}")
    project_id = args.project_id or root.name
    project_dir = args.project_dir or root / project_id
    output_path = args.output or root / "reports" / "sample_metrics.csv"
    file_count, row_count = compile_metrics(
        root, output_path, project_id, project_dir
    )
    print(f"Wrote {row_count} metric rows from {file_count} source CSVs to {output_path}")


if __name__ == "__main__":
    main()