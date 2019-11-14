StrainGE: Strain-level Genome Exploration
=========================================

StrainGE is a set of tools to analyse the within-species strain diversity in 
bacterial populations. It consists of two main components: 1) StrainGST: Strain
Genome Search tool, a tool to find close reference genomes for strains present
in a sample and 2) StrainGR: Strain Genome Recovery, a tool to perform
strain-aware variant calling at low coverages.

Dependencies
------------

* Python >= 3.6
* NumPy
* SciPy
* matplotlib
* scikit-bio
* pyvcf
* pysam
* h5py
* intervaltree

These dependencies will be automatically installed when creating a conda
environment from our `environment.yml` file described below.

Installation
------------

1. Install Anaconda or miniconda (if not already present on your system)
2. Clone the repository:

    `git clone https://github.com/broadinstitute/StrainGE`

3. Move into the directory:

    `cd StrainGE`

4. Create a new conda environment:

    `conda env create -f environment.yml`

5. Activate the environment:

    `source activate strainge`

6. Install StrainGE:

    `python setup.py install` or `python setup.py develop` for developer mode.

**TODO**: create our own pip/conda packages for even easier installation.

Running StrainGE
----------------

A typical workflow using StrainGE is as follows:

1. Build a database of known strains covering the diversity for a species of
   interest.
2. Run StrainGST to identify close reference genomes to strains present in
   a sample.
3. Run StrainGR to call SNPs and other variants, and to compare strains across
   samples.

Each step is divided into multiple subtasks, which can be ran using the
`strainge` command line program. The `strainge` program consists of multiple 
subcommands:

**Database creation**

* `strainge kmerize`
* `strainge kmersim`
* `strainge cluster`
* `strainge createdb`

**StrainGST: Strain Genome Search Tool**

* `strainge kmerize`
* `strainge search`

**StrainGR: Strain Genome Recovery**

* `strainge call`
* `strainge compare`
* `strainge dist`
* `strainge tree`

**Utilities**

* `strainge view`
* `strainge stats`
* `strainge plot`

This tutorial will explain how to run each step and how everything is 
connected. A graphical overview of the pipeline can be seen below.

![StrainGE pipeline](docs/img/strainge-overview.png)

### StrainGST database creation

#### 1. k-merize your reference sequences

This tutorial assumes you have activated the *strainge* conda environment.
Furthermore, we assume that you've downloaded your reference sequences as FASTA
files in the current working directory.

To k-merize all downloaded reference genomes, run the following bash command:

```bash
for f in *.fa; do straingst kmerize -k 23 -o $f.hdf5 $f; done;
```

For each FASTA file, there is now an accompanying HDF5 file containing the
k-mer data. With `-k` you can specify the k-mer size, which is by default 23.

#### 2. Compare the k-mer sets and cluster similar references

The goal of StrainGST is to identify close reference genomes to strains present
in a sample. These reference genomes are in turn used for variant calling and
sample comparisons. Here lies a trade-off: the reference genome should close
enough for accurate variant calling, but sample comparisons are more easy to
perform when the variant calling step is done using the same reference genome,
so you don't want to be too specific. The database of reference genomes should 
cover the diversity of the species of interest but not contain too many highly 
similar genomes. Therefore a clustering step is performed to reduce redundancy 
in the database.

We remove redundant reference genomes two ways:

1. Remove reference genomes that are a near perfect subset of another genome.
   An example of this an an *E. coli* strain used for synthetic biology
   applications that was basically a K-12 strain with many genes removed.
2. Cluster closely related genomes based on k-mer similarity and pick one 
   representative.

To do this, we need to compute the pairwise similarities between k-mer sets,
and a metric to identify whether a k-mer set is a subset of another. Both can
be obtained using `straingst kmersim`.

```bash
straingst kmersim --all-vs-all -t 4 -S jaccard -S subset *.hdf5 > similarities.tsv
```

This command produces as tab separated file, where each line contains
a pair of k-mer sets with their accompanying similarity scores. With the `-S`
flag we enable which scoring metrics to calculate, and in this case we enable
the *Jaccard* similarity and the *subset* score. The output file contains for 
each pair of k-mer sets the requested scores, sorted by the first scoring 
metric (in our case the jaccard similarity). With the parameter `-t` you 
specify the number of processes to spawn, to allow for parallel computation of 
these pairwise similarities.

We can now cluster our references using the `straingst cluster` command. 

```bash
straingst cluster -i similarities.tsv -d -C 0.99 -c 0.95 \
   --clusters-out clusters.tsv \
    *.hdf5 > references_to_keep.txt
```

The cluster command reads our previously created file `similarities.tsv` to
determine which references to keep. The first step is to discard any genome
where more than 99% of its kmers are present in another genome, as enabled by
`-d` and `-C 0.99`. Afterwards, we cluster similar genomes based on the
*Jaccard* similarity between k-mersets: if the Jaccard similarity between two
k-mer sets is higher than 0.95 (`-c 0.95`), those two genomes will be clustered 
together. For each cluster we pick one representative genome: the genome with 
the smallest mean distance to the other cluster members. Each genome to keep is
written to `references_to_keep.txt`. With the option `--clusters-out` we 
specify another file where we write the clustering results. Each line in this 
file specifies a cluster along with its entries, separated by a tab. The
genomes in the first column represent the cluster representatives. This option 
is optional, but can be useful for debugging purposes.

#### 3. Create pan-genome k-mer database

Using our list of references, we finally create a single database file which
will contain all k-mers of the given references. 

```bash
straingst createdb -f references_to_keep.txt -o pan-genome-db.hdf5
```

Now our database lives in the file `pan-genome-db.hdf5`, created from reference
sequences read from the file given by `-f`.

It is also possible to give the list of k-mer sets to include in the database
as positional arguments, like in the following example:

```bash
strainge createdb -o pan-genome-db.hdf5 ref1.hdf5 ref2.hdf5 ...
```

Combining the two methods described above works too.

### StrainGST: identify close reference genomes to strains in a sample

#### 1. k-merize your sample reads

With our database prepared, we can now analyse our samples, for example a 
metagenomic read data set from a patient. StrainGST uses k-mers to identify 
possible strains in your sample. Unsurprisingly, our first step is to kmerize 
our sample reads, for example if you have a FASTQ file named `patient1.fastq` 
with all the reads, then we generate its corresponding k-mer set as follows:

```bash
straingst kmerize -k 23 -o patient1.hdf5 patient1.fastq
```

Similar to the first step of the database creation section, this will generate 
a HDF5 file named `patient1.hdf5` with all k-mers and their corresponding 
counts. Make sure the value of `k` is the same as used in the 
database creation step.

You can specify multiple FASTQ files to the command above, which is useful if 
you have paired-end reads. Furthermore, it will also automatically decompress 
*gzipped* files. For example, if you have gzipped paired-end FASTQ files, then 
the following command also works:

```bash
straingst kmerize -k 23 -o patient1.hdf5 \
    patient1.1.fastq.gz patient1.2.fastq.gz
```

#### 2. Run StrainGST

We can now run `strainge search` with our database HDF5 and our sample HDF5:

```bash
straingst run -o results.tsv pan-genome-db.hdf5 patient1.hdf5
```

This will output a *tab separated values* (tsv) file, containing statistics 
about the sample k-mer set and a list of identified strains with accompanying 
statistics. For full description of all fields, please refer to the `straingst` 
reference documentation below.

### StrainGR: strain-aware variant calling and sample comparisons

#### 1. Inspect the StrainGST output and create a concatenated reference FASTA

Before variant calling it's required to align reads from a metagenomic sample 
to a reference. There are several options for what kind of reference to use:

1. A single reference genome from your database, for example the top strain
   identified by StrainGST. This will present difficulties
   if your sample has multiple strains, which will result in positions with
   evidence for multiple alleles. Furthermore, to perform variant calling
   comparisons across multiple samples, all samples will need to use the same
   reference, which may not be ideal.
2. Concatenate all reference genomes identified by StrainGST in a sample of 
   interest to a single file. This will introduce redundancy in the
   concatenated FASTA (because of shared sequence between strains), but the 
   aligner will be able to place reads with strain specific alleles to the 
   right location.
3. In many cases, you have multiple related samples (for example multiple
   samples from a patient at different time points). In many cases, you'll want
   to compare these samples with each other. In such cases it makes sense to
   collect all reference genomes identified by StrainGST across all these
   samples, and concatenate them to a single file.

StrainGR does not call variants in repetitive regions of the reference. By
concatenating genomes you will introduce repetitive regions if two strains both
have a collection of conserved genes for example. Reads from these regions will
placed ambiguously by the aligner. StrainGR detects this and will ignore those
locations. On the other hand, by concatenating reference genomes, reads with
strain-specific alleles will be placed at the correct location. In the end,
it's a trade off between how much of the genome will be ignored due to
repetitiveness and being able to place strain-specific reads correctly. 

It will depend on the project what option will be the best. The rest of the
tutorial assumes that a file `refs_concat.fasta` exists with the reference
genomes of interest.

#### 2. Align reads to the reference

StrainGR is built to be used with `bwa mem`, as it uses the supplied 
information on alternative alignment locations encoded in the `XA` SAM tag to 
deal with shared regions introduced by concatenating reference genomes.

The following command aligns the reads with `bwa mem` and outputs a sorted BAM
file:

```bash
bwa mem -t 2 refs_concat.fasta sample1.1.fq.gz sample1.2.fq.gz \
    | samtools sort -@ 2 -O BAM -o sample1.bam -

# Also create BAM index
samtools index sample1.bam
```

#### 3. Analyze read alignments to call variants

To call any variants in your sample run the StrainGR variant caller:

```bash
straingr call --hdf5-out sample1.hdf5 refs_concat.fasta sample1.bam
```

All variant calling data will be stored in the given HDF5 file
`sample1.hdf5`. When the variant caller is finished it will print a table with
summary statistics like coverage, SNP rate, gaps and more. You can also specify 
to output this table to a TSV file with the `-s` switch in the above command. 
There are more options for data output, it can output VCF files, BED tracks 
and more, see the CLI reference documentation below.

While the output to HDF5 is optional, it's strongly recommended to include this
output file. It will serve as input for sample comparisons, and all other 
output files (VCF, BED, summary stats) can be recreated at a later time using 
`straingr view`.


#### 4. Compare strains across samples

After running `straingr call` on a set of samples, it's possible to perform
strain-level comparisons. The following tools aid to analyze how closely
related strains across samples are:

* `straingr compare` - Nucleotide level comparisons of strains in different
   samples, outputs summary statistics.
* `straingr dist` - For all strains in samples close to a selected reference 
   genome, perform pairwise genetic distance estimation using Kimura's two
   parameter mode. Useful to create ordination plots.
* `straingr tree` - Using the distance matrix created by `straingr dist`,
   perform neighbour-joining to create an approximate phylogenetic tree.

CLI Reference Documentation
---------------------------

StrainGE has two main components, `straingst` and `straingr`, and each is
a separate program.

### Subcommand help

To view the documentation of each command line program, add the `-h` switch to
any subcommand. Examples:

```bash
straingst -h
straingst cluster -h
straingst run -h

straingr -h
straingr call -h
```

### Logging

StrainGE logs messages describing its current state to standard error. Its 
verboseness can be controlled using the `-v` switch, and this works for any of 
the commands described below. StrainGE has three verboseness levels. Please note 
that the `-v` switch needs to be set before any of the subcommands.

Examples:

```bash
straingst kmerize ...     # Only very global logging messages
straingst -v kmersim ...  # Get more internal logging messages
straingst -vv cluster ..  # Enable all debug messages
```
