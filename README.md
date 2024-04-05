# MAPseq Pipeline and Analysis [![DOI](https://zenodo.org/badge/741509600.svg)](https://zenodo.org/doi/10.5281/zenodo.10903736)

This repository provides an end-to-end open source pipeline for the
processing of MAPseq (**M**arker **A**ssisted **P**rogramming **seq**uencing) data.
* Current Version: v0.1

### Pre-Installation Requirements ###

We recommend creating a new conda environment for the installation of the tools required for the MAPseq pipeline and secondary analysis. For external software (not including Cell Ranger), we provide a [`requirements.yml`](requirements.yml) file for easy installation with conda. This can be run and initialized with:

```
conda env create -f requirements.yml
conda activate mapseq_env
```

#### External software required:

* Cell Ranger (v8.0.0) [(install here)](https://www.10xgenomics.com/support/software/cell-ranger/latest)
* Cell Ranger ATAC (v2.1.0) [(install here)](https://support.10xgenomics.com/single-cell-atac/software/pipelines/2.1/installation)
* kallisto (v0.46.1) [(install here)](https://pachterlab.github.io/kallisto/download)
* bustools (v0.39.3) [(install here)](https://bustools.github.io/download)
* cellbender (v0.3.0) [(install here)](https://cellbender.readthedocs.io/en/latest/installation/index.html)
* leidenalg (v0.10.0) [(install here)](https://leidenalg.readthedocs.io/en/stable/install.html)

### Installation ###
After installing necessary requirements, you may install the MAPseq pipeline from source using:

```
gh repo clone willisbillis/MAPseq
cd MAPseq
./install.sh
```

### Questions and Issues ###

If you have a question, error, or bug to report, please use the [issue page](https://github.com/willisbillis/MAPseq/issues).

Resources
---------
* [MAPseq tutorial](https://github.com/willisbillis/MAPseq/blob/main/docs/quickstart.md)

Citing
------
If you make use of this software for your work we would appreciate it if you would include the citation in any subsequent work. You may cite it using the dropdown menu in the top right of this page.

Known Issues
------------
* No currently known issues