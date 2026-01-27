# STEP 7: Testing & Final Documentation Plan

Status: In-progress | Target: End-to-end validation, performance checks, final docs

## Objectives
- Validate pipeline correctness on a small synthetic dataset (GEX + optional VDJ)
- Validate portal manifest integrity and DB mapping
- Capture performance benchmarks (GPU vs CPU) and resource footprints
- Finalize docs (testing section, benchmarking notes, troubleshooting updates)

## Test Scenarios
1) **Smoke (no VDJ)**
   - Inputs: tiny FASTQs (or use pre-aggregated H5AD checkpoints if available)
   - Params: `--harmony_enabled true`, `--include_vdj false`
   - Expected: QC → Harmony → Clustering outputs present; portal manifest lists clustered/harmony/annotated artifacts

2) **VDJ-enabled**
   - Inputs: include VDJ FASTQs; ensure `vdj_fastq_r1/r2` populated
   - Params: `--include_vdj true`
   - Expected: VDJ clonotype files + VDJ-integrated H5AD + VDJ metrics/summary in manifest

3) **CPU fallback**
   - Params: `--gpu_strategy off`, `--harmony_gpu false`, `--rapids_gpu false`
   - Expected: rapids fallbacks to scanpy; Harmony skips GPU; pipeline completes

4) **Portal ingest check**
   - After run, read `portal/<batch>/portal_manifest.json`
   - Validate md5/size for key files; confirm label/kind mapping
   - Optionally run a stub DB loader to insert rows into `files`/`metrics`

## Commands (examples)
```bash
# Run with minimal sample sheet (adjust paths)
nextflow run main.nf \
  --sample_sheet data/samples_tiny.csv \
  --ref_dir /path/to/ref \
  --outdir /tmp/mapseq_test \
  --genome GRCh38 \
  --gpu_strategy auto \
  --harmony_enabled true \
  --include_vdj false

# Re-run with VDJ enabled
nextflow run main.nf \
  --sample_sheet data/samples_tiny_vdj.csv \
  --include_vdj true \
  --outdir /tmp/mapseq_test_vdj
```

## Benchmarks to capture
- Wall-clock and max RSS for: CELLBENDER_QC, SCDBLFINDER_FILTER, HARMONY_INTEGRATION, RAPIDS_CLUSTERING
- GPU usage when enabled (nvidia-smi samples) vs CPU fallback timing
- Throughput: cells/sec for clustering step

## Validation checklist
- [ ] All expected artifacts exist for each batch (see PORTAL_INTEGRATION.md coverage)
- [ ] portal_manifest.json entries resolve to real files; md5/bytes match
- [ ] H5ADs load in scanpy without errors
- [ ] VDJ tables present when enabled; cluster summary not empty when clustering exists
- [ ] Metadata fields present: organism, donor_id, tissue, batch, sample_id

## Documentation tasks
- [ ] Add testing section to README/QUICK_REFERENCE with example commands
- [ ] Add benchmark results to DEPLOYMENT_SUMMARY.md (GPU vs CPU)
- [ ] Add troubleshooting entries for portal ingest (missing files/paths)

## Optional automation
- Add a `tests/` folder with a small manifest validator script that checks portal_manifest.json against filesystem and prints missing entries.
- Add a Makefile target: `make test-manifest OUTDIR=/tmp/mapseq_test BATCH=B001`
