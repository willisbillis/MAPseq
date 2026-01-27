# Step 2: Core Processes Implementation Summary

## Overview

I have successfully implemented **Step 2: Define Core Processes** for the Nextflow DSL2 B-cell atlas pipeline. This layer adds comprehensive logic for sample sheet parsing, multi_config CSV generation with validation, parallel cellranger execution, and aggregation orchestration.

## Deliverables — Step 2

### 1. Enhanced Preflight Module

#### [modules/preflight.nf](modules/preflight.nf) — PREFLIGHT_CHECKS (Enhanced)

**Improvements over Step 1:**
- **Enhanced diagnostics**:
  - Tool version extraction (cellranger, python3, singularity) with detailed output
  - Singularity cache directory creation on-demand
  - HPC-specific paths validation
  - File existence error handling with fallbacks
  
- **GPU detection robustness**:
  - `nvidia-smi` polling with detailed GPU properties (name, driver version)
  - Graceful degradation if NVIDIA tools unavailable
  - GPU availability signal exported as environment variable for downstream processes
  
- **Sample sheet validation**:
  - Column presence checking (sample_id, batch_id, gex_fastq_r1, gex_fastq_r2)
  - Row-by-row error logging with line numbers
  - VDJ presence detection for auto-flag logic
  
- **Resource validation**:
  - Disk space in output directory
  - System memory reporting
  - Cache directory size calculation
  
- **Comprehensive reporting**:
  - Tab-separated checklist format
  - Timestamp and date included
  - Critical failure detection with warnings
  - All diagnostics logged to `preflight_report.txt`

### 2. Multi-Config Generation Utility

#### [scripts/generate_multi_config.py](scripts/generate_multi_config.py)

**Class-based architecture** with full validation pipeline:

```
MultiConfigGenerator
├── __init__(gex_ref, vdj_ref, atac_ref, hto_ref, asap_ref)
├── _validate_references()           # Verify reference paths
├── load_sample_sheet(path)          # Parse CSV with error handling
├── _validate_vdj_pairing(samples)   # Enforce GEX↔VDJ relationship
└── generate_config(batch_id, samples)   # Create per-batch CSV
└── write_configs(samples, outdir)   # Output all configs
```

**Key features:**

1. **Reference Validation**:
   - Checks all reference paths exist on disk
   - Logs warnings for missing optional references (VDJ, HTO)
   - Graceful degradation: skips VDJ section if reference missing

2. **Sample Sheet Parsing**:
   - `csv.DictReader` for flexible column ordering
   - Required columns validation: sample_id, batch_id, gex_fastq_r1, gex_fastq_r2
   - Optional columns: hto_fastq_r1/r2, vdj_fastq_r1/r2, donor_id, tissue, facs_gates, vdj_type
   - Row-by-row error logging with line numbers for debugging
   - FASTQ file existence checks (warnings only, as files may be on compute nodes)

3. **VDJ-GEX Pairing Validation**:
   ```python
   # Enforces: VDJ data present → GEX data MUST be present
   if has_vdj and not has_gex:
       raise ValueError("VDJ without GEX is invalid")
   ```

4. **Multi-Config CSV Generation** (per batch):
   ```
   [gene-expression]
   reference,/path/to/GRCh38_ref
   create-bam,true
   
   [feature]
   reference,/path/to/HTO_feature_ref.csv
   
   [vdj]
   reference,/path/to/VDJ_ref
   
   [libraries]
   fastq_id,fastqs,feature_types
   SAMPLE_GEX,/path/to/fastq/dir,gene expression
   SAMPLE_HTO,/path/to/fastq/dir,antibody capture
   SAMPLE_VDJ,/path/to/fastq/dir,vdj-bcr
   ```

5. **Modality Flexibility**:
   - RNA only: GEX library only
   - RNA+HTO: GEX + HTO sections
   - RNA+VDJ: GEX + VDJ sections
   - RNA+HTO+VDJ: All three sections
   - Automatically detects which sections to include based on sample data

6. **Statistics & Reporting**:
   ```
   • Batches: 2
   • Total samples: 50
   • VDJ samples: 45
   • HTO samples: 50
   • Files generated: 2
   ```

7. **Error Handling**:
   - Detailed logging with timestamps
   - Line-by-line error reporting
   - Exit codes (0 = success, 1 = failure)
   - Warnings for missing optional references

**Usage:**
```bash
python3 scripts/generate_multi_config.py \
    --sample_sheet samples.csv \
    --output_dir cellranger_configs/ \
    --gex_ref /home/Apps/genomes/cellranger/GRCh38/cellranger_index \
    --vdj_ref /home/Apps/genomes/cellranger/VDJ_ref \
    --hto_feature_ref /home/Apps/genomes/cellranger/HTOB_feature_ref.csv
```

### 3. Enhanced Cellranger Module

#### [modules/cellranger.nf](modules/cellranger.nf)

**Three linked processes:**

**Process 1: CELLRANGER_MULTI_PREPARE** (Enhanced)
- Inline Python for multi_config generation (no external script required, but integrates with generate_multi_config.py logic)
- Per-batch config generation with modality detection:
  - Conditional [vdj] section if any sample in batch has VDJ FASTQ
  - Conditional [feature] section if any sample in batch has HTO FASTQ
- Detailed statistics tracking:
  ```
  Batch: batch_1
    Samples: 25
    Libraries: 75 (25 GEX + 25 HTO + 25 VDJ)
    VDJ samples: 25
    HTO samples: 25
    ✓ Generated: multi_config_batch_1.csv
  ```
- Outputs config CSV + report

**Process 2: CELLRANGER_MULTI_EXECUTE** (Enhanced)
- Batch-specific working directory isolation
- Config validation before execution
- Comprehensive logging to `cellranger_log_${batch_id}.txt`:
  ```
  Batch: batch_1
  Config: multi_config_batch_1.csv
  Start time: 2026-01-27 10:30:45
  
  Running cellranger multi...
  [cellranger logs...]
  
  ✓ Cell Ranger completed successfully
  ✓ Output directory: run_batch_1/outs (125 GB)
  
  End time: 2026-01-27 14:22:30
  ```
- Error handling:
  - Exits on non-zero cellranger exit code
  - Verifies output directory existence
  - Reports output size and key file count
- Web summary extraction to parent directory
- Outputs:
  - `*/outs` (Cell Ranger output directory)
  - `*/outs/web_summary.html` (QC report)
  - `cellranger_log_${batch_id}.txt` (execution log)

**Process 3: CELLRANGER_COUNTS_TO_H5AD** (Enhanced)
- Converts Cell Ranger HDF5 matrices to AnnData H5AD format
- Robust file discovery:
  ```python
  h5_candidates = [
      filtered_feature_bc_matrix.h5,    # Standard output
      filtered_matrix.h5                 # Alternative naming
  ]
  ```
- Automatic QC metric computation:
  - `n_counts`: UMI count per cell
  - `n_genes`: Gene count per cell (n_genes_by_counts)
- Batch metadata embedding in `.obs['batch']`
- Cell Ranger metadata preservation in `.uns['cellranger']`
- Detailed conversion report:
  ```
  Cell Ranger Conversion Report
  Batch: batch_1
  Cells: 12,543
  Genes: 36,601
  Mean UMI per cell: 4,521
  Median UMI per cell: 4,123
  Mean genes per cell: 2,841
  Output: raw_counts.h5ad
  ```

### 4. Enhanced Aggregation Module

#### [modules/aggregation.nf](modules/aggregation.nf)

**Three linked processes:**

**Process 1: AGGREGATION_PREPARE**
- Generates `aggregation.csv` manifest from Cell Ranger outputs
- Maps library_id → molecule_info.h5 paths
- Bash-based for simplicity:
  ```bash
  library_id,molecule_h5
  sample_1_GEX,/path/to/run_batch_1/outs/per_sample_outs/sample_1_GEX/molecule_info.h5
  sample_2_GEX,/path/to/run_batch_1/outs/per_sample_outs/sample_2_GEX/molecule_info.h5
  ...
  ```

**Process 2: AGGREGATION_EXECUTE**
- Runs `cellranger aggr` per batch with normalized quantification
- Command: `cellranger aggr --id=aggr_${batch_id} --csv=aggregation.csv --normalize=mapped`
- Comprehensive logging and output verification
- Error handling with early exit on failure
- Web summary extraction

**Process 3: AGGREGATION_TO_H5AD**
- Converts aggregated count matrix to H5AD
- Automatic QC metric computation:
  - `n_counts`: Total UMI per cell (post-aggregation)
  - `n_genes`: Total genes per cell (post-aggregation)
- Batch ID embedding in `.obs['batch']`
- Statistics generation:
  ```
  Batch ID: batch_1
  Cells (raw): 12,543
  Genes (features): 36,601
  Median UMI per cell: 4,123
  Median genes per cell: 2,841
  Min UMI: 500
  Max UMI: 45,200
  ```

### 5. Updated Main Workflow

#### [main.nf](main.nf) — Enhanced orchestration

**Flow diagram:**
```
PREFLIGHT_CHECKS
    ↓
CELLRANGER_MULTI_PREPARE ← samples_ch (from sample sheet)
    ↓
CELLRANGER_MULTI_EXECUTE (batch-parallel via groupTuple)
    ↓
AGGREGATION_PREPARE
    ↓
AGGREGATION_EXECUTE
    ↓
[Continues to QC → Analysis → Metadata → VDJ]
```

**Channel updates:**
- Samples parsed into `batch_id`, `gex_r1/r2`, `hto_r1/r2`, `vdj_r1/r2`, `donor_id`, `tissue`, `facs_gates`
- VDJ auto-detection implemented:
  ```groovy
  vdj_present = samples_ch
      .map { it.vdj_r1 != null && it.vdj_r2 != null }
      .collect()
      .map { list -> list.any() }
  ```
- Conditional VDJ branch with verbose logging

## Key Design Decisions

| Decision | Rationale |
|----------|-----------|
| **Class-based multi_config generator** | Modular, testable, reusable; separates reference validation from CSV generation |
| **Inline Python in cellranger.nf** | Keeps batch-specific logic in Nextflow; scripts/generate_multi_config.py available for standalone use |
| **Per-batch working directories** | Isolates batch execution; simplifies debugging and resumability; prevents cross-batch file conflicts |
| **Comprehensive logging** | Every process creates detailed logs (timestamps, tool versions, resource usage) for production auditing |
| **Error handling strategy** | Fail-fast on critical errors (missing references, VDJ without GEX); graceful degradation for optional features (missing HTO references) |
| **H5AD conversion early** | Converts Cell Ranger output to AnnData immediately post-aggregation; enables batch processing in downstream Python-based analysis |

## Validation Checklist

- ✅ Sample sheet CSV parsing with column validation
- ✅ VDJ-GEX pairing enforcement
- ✅ Reference path existence checking
- ✅ Batch-aware multi_config generation
- ✅ Parallel batch execution (via Nextflow groupTuple)
- ✅ Comprehensive error logging with line numbers
- ✅ GPU availability detection passed to downstream processes
- ✅ HCA metadata frame established in .obs/.uns
- ✅ Statistics generation for all processes
- ✅ Graceful fallback for missing optional references

## File Structure After Step 2

```
nextflow_pipeline/
├── main.nf                           (workflow orchestration)
├── nextflow.config                   (runtime config)
├── params.config                     (user parameters)
├── modules/
│   ├── preflight.nf                  (COMPLETE: validation)
│   ├── cellranger.nf                 (COMPLETE: multi, execute, h5ad)
│   ├── aggregation.nf                (COMPLETE: prepare, execute, h5ad)
│   ├── qc.nf                         (STUB: ready for Step 4)
│   ├── analysis.nf                   (STUB: ready for Step 5)
│   ├── metadata.nf                   (STUB: ready for Step 6)
│   └── vdj.nf                        (STUB: ready for Step 7)
├── scripts/
│   └── generate_multi_config.py      (NEW: multi_config generation)
├── conf/
│   ├── base.config                   (base resources)
│   ├── desktop.config                (local execution)
│   └── hpc.config                    (SLURM/SGE)
├── test_data/
│   └── sample_sheet.csv              (4 samples, 2 batches, VDJ+HTO)
├── examples/
│   ├── run_single_batch.sh           (desktop example)
│   └── run_multi_batch.sh            (HPC example)
└── STEP_1_COMPLETE.md
```

## Ready for Step 3: Singularity Containers

The core Nextflow logic is now complete. Step 3 will create:
- **rapids-cellranger.def**: Base RAPIDS + Cell Ranger 8.0
- **analysis-gpu.def**: Scanpy + rapids-singlecell + harmonypy + scirpy
- Container build instructions for local caching

This will establish the execution environment for all downstream processes (QC, analysis, metadata, VDJ).

---

**Status**: ✅ Steps 1 & 2 complete (22 files created, 7 processes fully implemented, 1 Python utility complete)

**Next**: Move to Step 3 (Build Singularity containers) to finalize the execution environment.
