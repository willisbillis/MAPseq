# MAPseq Pipeline and Analysis [![DOI](https://zenodo.org/badge/741509600.svg)](https://zenodo.org/doi/10.5281/zenodo.10903736)

![logo](docs/logo.png)

This repository provides an end-to-end open source pipeline for the
processing of MAPseq (**M**arker **A**ssisted **P**rogramming **seq**uencing) data.

## Pre-Installation Requirements

Supported Operating Systems:
 - Ubuntu 22.04


We recommend creating a new conda environment for the installation of the tools required for the MAPseq pipeline and secondary analysis. For external software (not including Cell Ranger), we provide an `environment.yml` file for easy installation with conda. This can be run and initialized with:

# Clone the repository
git clone https://github.com/willisbillis/MAPseq.git
# Or, if you have GitHub CLI installed:
# gh repo clone willisbillis/MAPseq

cd MAPseq
conda env create --name mapseq_env --file environment_<your_machine_os>.yml
conda activate mapseq_env
```

### External software required

* Cell Ranger (v8.0.0) [(install here)](https://www.10xgenomics.com/support/software/cell-ranger/latest)
* Cell Ranger ATAC (v2.1.0) [(install here)](https://support.10xgenomics.com/single-cell-atac/software/pipelines/2.1/installation)
* kallisto (<v1.0) [(install here)](https://pachterlab.github.io/kallisto/download)
* bustools (<v1.0) [(install here)](https://bustools.github.io/download)
* cellbender (v0.3.0) [(install here)](https://cellbender.readthedocs.io/en/latest/installation/index.html)
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

Citing

---------

If you make use of this software for your work we would appreciate it if you would include the citation in any subsequent work. You may cite it using the dropdown menu in the top right of this page.

Known Issues
------------
* No currently known issues
