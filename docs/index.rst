.. StrainGE documentation master file, created by
   sphinx-quickstart on Thu Jan 14 11:54:13 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

========
StrainGE
========
A toolkit to track and characterize low-abundance strains using metagenomic data
--------------------------------------------------------------------------------



StrainGE is a set of tools to analyse conspecific strain diversity in 
bacterial populations. It consists of two main components: 1) Strain
Genome Search tool (StrainGST), a tool to find close reference genomes to strain(s) present
in a sample and 2) Strain Genome Recovery (StrainGR), a tool to perform
strain-aware variant calling at low coverages, which in turn can be used to track strains across samples.

.. image:: img/strainge-overview-full.png
   :width: 500px

Installation
------------

StrainGE requires Python >= 3.7 and depends on the following packages:

* NumPy
* SciPy
* matplotlib
* scikit-bio
* pyvcf
* pysam
* h5py
* intervaltree

These packages will be automatically installed when installing through pip.

Install through `pip`
=====================

.. code-block:: bash

   pip install strainge

Install manually from github
============================

1. Clone the repository
   
   .. code-block:: bash

      git clone https://github.com/broadinstitute/StrainGE

2. Install StrainGE

   .. code-block:: bash

      cd StrainGE
      python setup.py install

Usage
-----

.. toctree::
   :maxdepth: 1

   create-database
   straingst
   straingr
   compare-strains


Citation
--------

If you use StrainGE in your project, please consider citing our publication: TODO

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
