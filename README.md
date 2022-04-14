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

**Warning**: NumPy already has to be installed otherwise the above command will fail.
You'll have to make sure all tools like `bwa`, `samtools` and `mummer` are installed as well.

### Conda

1. Install Anaconda or miniconda (if not already present on your system)
2. Create a new environment:

    `conda create -n strainge python=3`
3. Activate the environment:

    `source activate strainge`
4. Enable `bioconda` and `conda-forge` channels:

   ```
   conda config --add channels bioconda
   conda config --add channels conda-forge
   ```
5. Install StrainGE:
   
   `conda install strainge`

Optional tip: also consider installing [mamba](https://github.com/mamba-org/mamba) before installing StrainGE for much 
faster conda operations.


Documentation
-------------

The documentation can be read on [readthedocs](https://strainge.readthedocs.io).

Citation
--------

Dijk, Lucas R. van, Bruce J. Walker, Timothy J. Straub, Colin J. Worby, Alexandra Grote, Henry L. Schreiber, Christine Anyansi, et al. 2022. “StrainGE: A Toolkit to Track and Characterize Low-Abundance Strains in Complex Microbial Communities.” Genome Biology 23 (1): 74. https://doi.org/10.1186/s13059-022-02630-0.

