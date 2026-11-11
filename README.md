# MAPseq Pipeline and Analysis [![DOI](https://zenodo.org/badge/741509600.svg)](https://zenodo.org/doi/10.5281/zenodo.10903736)

![logo](assets/logo.png)

This repository provides an end-to-end open source pipeline for the
processing of MAPseq (**M**arker **A**ssisted **P**rogramming **seq**uencing) data.

## Pre-Installation Requirements

Supported Operating Systems:
 - Ubuntu 24.04


We recommend creating a new conda environment for the installation of the tools required for the MAPseq pipeline and secondary analysis. For external software (not including Cell Ranger), we provide an `environment.yml` file for easy installation with conda. This can be run and initialized with:

```bash
# Clone the repository
git clone https://github.com/willisbillis/MAPseq.git
# or, if you have GitHub CLI installed:
gh repo clone willisbillis/MAPseq

cd MAPseq
conda env create --name mapseq_env --file environment_<your_machine_os>.yml
conda activate mapseq_env
```

### External software required

* Cell Ranger (v10.0.0) [(install here)](https://www.10xgenomics.com/support/software/cell-ranger/latest)
* Cell Ranger ATAC (v2.1.0) [(install here)](https://www.10xgenomics.com/support/software/cell-ranger-atac/downloads/installation)
* kb (kallisto | bustools) (<v1.0) [(install here)](https://kallisto.readthedocs.io/en/latest/quick_start/installation.html)
* leidenalg (<v1.0) [(install here)](https://leidenalg.readthedocs.io/en/stable/install.html)

### Installation

After installing necessary requirements, you may install the MAPseq pipeline from source using:

```bash
gh repo clone willisbillis/MAPseq
cd MAPseq
./install.sh # sudo password will be prompted during this script
```

It is necessary to reload the shell for changes to take effect after installation. We recommend doing this by closing the current session and opening a new one.

### Remove an Installation

To remove changes made by the installation, you may navigate to the MAPseq repository and revert the changes with the following:

```bash
cd $MAPSEQ_REPO_PATH
./uninstall.sh # sudo password will be prompted during this script
```

For the changes to take effect, the shell must be reloaded. We recommend doing this by closing the current session and opening a new one.

### Questions and Issues

If you have a question, error, or bug to report, please use the [issue page](https://github.com/willisbillis/MAPseq/issues).

Resources
---------

* [MAPseq tutorial](docs/quickstart.md)

## Pipeline Metrics

Mapping keeps the Cell Ranger HTML summaries and also copies the native metrics
CSV into each modality's `pipeline/*/reports/` directory. The aggregation
driver compiles available mapping, demultiplexing, and QC CSVs automatically at
the end of an aggregation run. Since QC scripts run separately, rerun the
compiler after QC to refresh the ledger with those results:

```bash
python3 /path/to/aggregation_run/phase2_scripts/compile_pipeline_metrics.py /path/to/project_parent
```

The output is `reports/sample_metrics.csv`, a long-form table with one row per
metric and columns for project, run, modality, stage, sample ID, sample scope,
value, and source file. `sample_scope` distinguishes Cell Ranger samples,
demultiplexing libraries, and QC sample IDs. Existing stage-specific CSVs
remain available alongside this combined ledger; the HTML reports are retained
for review.

## Generic Phase-2 QC

New aggregation projects also include configurable RNA and ATAC phase-2
templates in `phase2_scripts/`. Set the project paths, `MAPSEQ_REPO_PATH`, sample
ID columns, and QC thresholds at the top of each script, inspect the generated
QC plots, then run the corresponding script with `Rscript`. The helper is loaded
from the local MAPseq repository, so the templates can be moved without moving
the helper or relying on their location. They write a timestamped backup of
`reports/sample_metrics.csv` before appending modality-specific QC metrics, plus
Spearman correlation statistics and plots for filtered sample cell count versus
mean sequencing depth. Re-running replaces that template's prior QC rows instead
of duplicating them.

The templates follow the current [single-cell best-practices guidance on QC](https://www.sc-best-practices.org/preprocessing-visualization/quality-control/),
[feature selection](https://www.sc-best-practices.org/preprocessing-visualization/feature-selection/),
[dimensionality reduction](https://www.sc-best-practices.org/preprocessing-visualization/dimensionality-reduction/),
and [ATAC QC and embedding](https://www.sc-best-practices.org/chromatin-accessibility/dimensionality-reduction-clustering/).
Thresholds are dataset-specific, so the scripts produce joint QC plots for
review rather than treating defaults as universal. The ATAC template keeps
Signac TF-IDF/LSI to match the MAPseq toolchain and makes fragment-depth effects
visible on the embedding; the guide discusses SnapATAC2 spectral embedding as
an alternative that would require a separate toolchain.

Citing
---------

If you make use of this software for your work we would appreciate it if you would include the citation in any subsequent work. You may cite it using the dropdown menu in the top right of this page.
