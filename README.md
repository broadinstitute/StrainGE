StrainGE: Strain-level Genome Exploration
=========================================

StrainGE is a set of tools to analyse the within-species strain diversity in 
bacterial populations. It consists of two main components: 1) StrainGST: Strain
Genome Search tool, a tool to find close reference genomes for strains present
in a sample and 2) StrainGR: Strain Genome Recovery, a tool to perform
strain-aware variant calling at low coverages.

[![Documentation Status](https://readthedocs.org/projects/strainge/badge/?version=latest)](https://strainge.readthedocs.io/en/latest/?badge=latest)
[![Python Package Index](https://img.shields.io/pypi/v/strainge)](https://pypi.org/project/strainge)
[![DOI](https://zenodo.org/badge/71356920.svg)](https://zenodo.org/badge/latestdoi/71356920)



Dependencies
------------

### Python packages

* Python >= 3.7
* NumPy
* SciPy
* matplotlib
* scikit-bio >= 0.5
* scikit-learn >= 0.24
* pyvcf
* pysam
* h5py
* intervaltree

### Bioinformatics tools

* bwa
* samtools
* mummer


Installation
------------

### Python Package Index

`pip install strainge`

You'll have to make sure all tools like `bwa`, `samtools` and `mummer` are installed as well.

### Conda

1. Install Anaconda or miniconda (if not already present on your system)
2. Clone the repository:

    `git clone https://github.com/broadinstitute/StrainGE`

3. Move into the directory:

    `cd StrainGE`

4. Create a new conda environment:

    `conda env create -f environment.yml`

5. Activate the environment:

    `source activate strainge`


Documentation
-------------

The documentation can be read on [readthedocs](https://strainge.readthedocs.io).

Citation
--------

Dijk LR van, Walker BJ, Straub TJ, Worby C, Grote A, Schreiber HL, et al. StrainGE: A toolkit to track and characterize low-abundance strains in complex microbial communities. bioRxiv. 2021 Feb 14;2021.02.14.431013. 
