StrainGE: Strain-level Genome Exploration
=========================================

StrainGE is a set of tools to analyse the within-species strain diversity in 
bacterial populations. It consists of two main components: 1) StrainGST: Strain
Genome Search tool, a tool to find close reference genomes for strains present
in a sample and 2) StrainGR: Strain Genome Recovery, a tool to perform
strain-aware variant calling at low coverages.

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
3. Run StrainGR to call SNPs and other variants

Each step is divided into multiple subtasks, which can be ran using the
`strainge` command line program. The `strainge` program consists of multiple 
subcommands:

**Database creation**

* `strainge kmerize`
* `strainge compare`
* `strainge cluster`
* `strainge createdb`

**StrainGST: Strain Genome Search Tool**

* `strainge kmerize`
* `strainge search`

**StrainGR: Strain Genome Recovery**

* `strainge call`
* `strainge sample-compare`

**Utilities**

* `strainge view`
* `strainge stats`
* `strainge plot`

This tutorial will explain how to run each step and how everything is 
connected. A graphical overview of the pipeline can be seen below.

![StrainGE pipeline](docs/img/strainge-overview.png)

### Database creation

#### 1. k-merize your reference sequences

This tutorial assumes you have activated the *strainge* conda environment.
Furthermore, we assume that you've downloaded your reference sequences as FASTA
files in the current working directory.

To k-merize all downloaded reference genomes, run the following bash command:

```bash
for f in *.fa; do strainge kmerize -k 23 -o $f.hdf5 $f; done;
```

For each FASTA file, there is now an accompanying HDF5 file containing the
k-mer data. With `-k` you can specify the k-mer size, which is by default 23.

#### 2. Compare the k-mer sets and cluster similar references

The goal of StrainGST is to identify close reference genomes to strains present
in a sample. These reference genomes are in turn used for variant calling and
sample comparisons. Here lies a trade-off: the reference genome should close
enough for accurate variant calling, but sample comparisons are more easy to
perform when the variant calling step is done using the same reference genome. 
The database of reference genomes should cover the diversity of the
species of interest but not contain too many highly similar genomes. Therefore
a clustering step is performed to reduce redundancy in the database.

First, we need to compute the pairwise similarities between k-mer sets. This
can be done using `strainge compare`:

```bash
strainge compare --all-vs-all -t 4 *.hdf5 > similarities.tsv
```

This command produces as tab separated file, where each line contains
a pair of k-mer sets with their accompanying similarity scores. The file is
sorted such that highly similar pairs come first. With the parameter `-t` you 
specify the number of processes to spawn, to allow for parallel computation of 
these pairwise similarities.

We can now cluster our references using the `strainge cluster` command. 

```bash
strainge cluster -i similarities.tsv -c 0.95 --clusters-out clusters.tsv \
    *.hdf5 > references_to_keep.txt
```

The cluster command reads our previously created file `similarities.tsv` and
figures out which references are more than 95% similar to each other, as
specified with `-c`. Within each cluster, the reference with the least amount
of scaffolds is kept, and its filename will be written to standard output. We
store the references we want to keep in the file `references_to_keep.txt`. With
the option `--clusters-out` we specify another file where we write the
clustering results. Each line in this file specifies a cluster along with each
entries, separated by a tab. The entries within a cluster are also sorted by
number of scaffolds in the FASTA file, which means that the first entry on this
line is the reference with the least amount of scaffolds, and in other words is
the reference that will be included in the database. The rest of the entries 
are ignored. This option is optional, but can be useful for debugging purposes.

#### 3. Create pan-genome k-mer database

Using our list of references, we finally create a single database file which
will contain all k-mers of the given references. 

```bash
strainge createdb -f references_to_keep.txt -o pan-genome-db.hdf5
```

Now our database lives in the file `pan-genome-db.hdf5`, created from reference
sequences read from the file given by `-f`.

It is also possible to give the list of k-mer sets to include in the database
as positional arguments, like in the following example:

```bash
strainge createdb -o pan-genome-db.hdf5 ref1.hdf5 ref2.hdf5 ...
```

Combining the two methods described above works too.

#### Appendix A. Combining above steps

Using bash we can pipe a lot of data directly to the next step. So to reduce
disk writes, we can combine almost all steps described above to one single bash
command:

```bash
strainge compare --all-vs-all -t 4 *.hdf5 \
    | strainge cluster -i - -c 0.95 --clusters-out clusters.tsv *.hdf5 \
    | strainge createdb -F -f - -o pan-genome-db.hdf5
```

By using `-` we signify that `strainge` should read from the standard input
instead of a file from disk.

### StrainGST: identify close reference genomes to strains in a sample

#### 1. k-merize your sample reads

With our database prepared, we can now analyse our samples, for example a 
metagenomic read data set from a patient. StrainGST uses k-mers to identify 
possible strains in your sample. Unsurprisingly, our first step is to kmerize 
our sample reads, for example if you have a FASTQ file named `patient1.fastq` 
with all the reads, then we generate its corresponding k-mer set as follows:

```bash
strainge kmerize -k 23 -o patient1.hdf5 patient1.fastq
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
strainge kmerize -k 23 -o patient1.hdf5 \
    patient1.1.fastq.gz patient1.2.fastq.gz
```

#### 2. Run StrainGST

We can now run `strainge search` with our database HDF5 and our sample HDF5:

```bash
strainge search -o results.tsv pan-genome-db.hdf5 patient1.hdf5
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
   interest to a single file. This will introduce redundancy in the reference 
   FASTA (because of shared sequence between strains), but the aligner will 
   be able to place reads with strain specific alleles to the right location.
3. In many cases, you have multiple related samples (for example multiple
   samples from a patient at different time points). In many cases, you'll want
   to compare these samples with each other. In such cases it makes sense to
   collect all reference genomes identified by StrainGST across all these
   samples, and concatenate them to a single file.

In the end, it's a trade-off between alignment accuracy of reads unique to
a strain, alignment accuracy of reads from shared regions, and
how easy it is to compare samples (which requires that reads are aligned to the
same FASTA file).

It will depend on the project what option will be the best. The rest of the
tutorial assumes that a file `refs_concat.fasta` exists with the reference
genomes of interest.

#### 2. Align reads to the reference

In theory, StrainGR could work with any aligner. In practice, however, it's
recommended to use `bwa mem`, as it uses the supplied information on alternative
alignment locations encoded in the `XA` SAM tag to deal with shared regions
introduced by concatenating reference genomes.

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
strainge call --hdf5-out sample1.hdf5 refs_concat.fasta sample1.bam
```

All variant calling data will be stored in the given HDF5 file
`sample1.hdf5`. When the variant caller is finished it will print a table with
summary statistics like coverage, SNP rate, gaps and more. You can also specify 
to output this table to a TSV file with the `-s` switch in the above command. 
There are more options for data output, it can output VCF files, BED tracks 
and more, see the CLI reference documentation below.

While the output to HDF5 is optional, it's strongly recommended to include this
output file. It will serve as input for sample comparisons, and all other output
files (VCF, BED, summary stats) can be recreated at a later time using `strainge
view`.

CLI Reference Documentation
---------------------------

The main CLI entry point is the program called `strainge`. This tool is divided
in to multiple subcommands described below.

StrainGE logs messages describing its current state to standard error. Its 
verboseness can be controlled using the `-v` switch, and this works for any of 
the commands described below. StrainGE has three verboseness levels. Please note 
that the `-v` switch needs to be set before any of the subcommands.

Examples:

```bash
strainge kmerize ...     # Only very global logging messages
strainge -v compare ...  # Get more internal logging messages
strainge -vv cluster ..  # Enable all debug messages
```

### Database management

#### `strainge kmerize` - kmerize references and samples

    usage: strainge kmerize [-h] [-k K] [-o OUTPUT] [-f] [-s SKETCH_FRACTION] [-F]
                            [-l LIMIT] [-p PRUNE]
                            sequences [sequences ...]

    K-merize a given reference sequence or a sample read dataset.

    positional arguments:
      sequences             Input sequence files (fasta or fastq by default;
                            optionally compressed with gz or bz2)

    optional arguments:
      -h, --help            show this help message and exit
      -k K, --k K           K-mer size (default 23)
      -o OUTPUT, --output OUTPUT
                            Filename of the output HDF5.
      -f, --fingerprint     Compute and save min-hash fingerprint (sketch).
      -s SKETCH_FRACTION, --sketch-fraction SKETCH_FRACTION
                            Fraction of k-mers to keep for a minhash sketch.
                            Default: 0.01
      -F, --filter          Filter output kmers based on kmer spectrum (to prune
                            sequencing errors)
      -l LIMIT, --limit LIMIT
                            Only process about this many kmers (can have suffix of
                            M or G)
      -p PRUNE, --prune PRUNE
                            Prune singletons after accumulating this (can have
                            suffix of M or G)

#### `strainge compare` - compare k-mer sets 

    usage: strainge compare [-h] (-a | -s FILE) [-f]
                            [-S {jaccard,minsize,meansize,maxsize,reference}] [-F]
                            [-t THREADS] [-o FILE]
                            strains [strains ...]

    Compare k-mer sets with each other. Both all-vs-all and one-vs-all is
    supported.

    positional arguments:
      strains               Filenames of k-mer set HDF5 files.

    optional arguments:
      -h, --help            show this help message and exit
      -a, --all-vs-all      Perform all-vs-all comparisons for the given k-mer
                            sets. Either --all-vs-all is required or --sample.
      -s FILE, --sample FILE
                            Perform one-vs-all comparisons with the given filename
                            as sample. Either --all-vs-all is required or
                            --sample.
      -f, --fingerprint     Use min-hash fingerprint instead of full k-mer set.
      -S {jaccard,minsize,meansize,maxsize,reference}, --scoring {jaccard,minsize,meansize,maxsize,reference}
                            The scoring metric to use (default: jaccard).
      -F, --fraction        Output numerator and denominator separately instead of
                            evaluating the division.
      -t THREADS, --threads THREADS
                            Use multiple processes the compute the similarity
                            scores (default 1).
      -o FILE, --output FILE
                            File to write the results (default: standard output).

#### `strainge cluster` - cluster similar k-mersets

    usage: strainge cluster [-h] [-c CUTOFF] [-i FILE] [-o FILE]
                            [--clusters-out FILE]
                            kmerset [kmerset ...]

    Group k-mer sets that are very similar to each other together.

    positional arguments:
      kmerset               The list of HDF5 filenames of k-mer sets to cluster.

    optional arguments:
      -h, --help            show this help message and exit
      -c CUTOFF, --cutoff CUTOFF
                            Minimum similarity between two sets to group them
                            together.
      -i FILE, --similarity-scores FILE
                            The file with the similarity scores between sets (the
                            output of 'strainge compare --all-vs-all'). Defaults
                            to standard input.
      -o FILE, --output FILE
                            The file where the list of sets to keep after
                            clustering gets written. Defaults to standard output.
      --clusters-out FILE   Output an optional tab separated file with all
                            clusters and their entries.

#### `strainge createdb` - create a pan-genome database from multiple k-mer sets

    usage: strainge createdb [-h] [-o OUTPUT] [-F] [-f FILE]
                             [kmerset [kmerset ...]]

    Create pan-genome database in HDF5 format from a list of k-merized
    strains.

    positional arguments:
      kmerset               The HDF5 filenames of the kmerized reference strains.

    optional arguments:
      -h, --help            show this help message and exit
      -o OUTPUT, --output OUTPUT
                            Pan-genome database output HDF5 file.
      -F, --fingerprint     Create fingerprint for the pan-genome database.
      -f FILE, --from-file FILE
                            Read list of HDF5 filenames to include in the database
                            from a given file (use '-' to denote standard input).
                            This is in addition to any k-merset given as
                            positional argument.

### StrainGST: strain genome search tool

#### `strainge search` - identify close reference genomes to strains in a sample

    usage: strainge search [-h] [-o OUTPUT] [-i ITERATIONS] [-t TOP] [-f]
                           [-F MINFRAC] [-s SCORE] [-e EVENNESS]
                           pan sample

    StrainGST: strain genome search tool. Identify close reference genomes
    to strains present in a sample.

    positional arguments:
      pan                   hdf5 file containing pan genome kmer set
      sample                Search for strains in this sample

    optional arguments:
      -h, --help            show this help message and exit
      -o OUTPUT, --output OUTPUT
                            output text file (default: standard out)
      -i ITERATIONS, --iterations ITERATIONS
                            max strains to look for (default: 5)
      -t TOP, --top TOP     How many best matches to print per iteration (default:
                            1)
      -f, --fingerprint     Using fingerprint rather than whole kmer set
      -F MINFRAC, --minfrac MINFRAC
                            minimum fraction of original kmers in strain (default:
                            0.01)
      -s SCORE, --score SCORE
                            minimum score (default: 0.01)
      -e EVENNESS, --evenness EVENNESS
                            minimum evenness (default: 0.60)

### StrainGR: strain genome recovery

#### `strainge call` - strain-aware variant caller

    usage: strainge call [-h] [-Q MIN_QUAL] [-P MIN_PILEUP_QUAL]
                         [-F MIN_QUAL_FRAC] [-M MIN_MAPPING_QUAL] [-G MIN_GAP]
                         [-s FILE] [--hdf5-out FILE] [-V FILE] [--verbose-vcf]
                         [--track-callable FILE] [--track-poor-mq FILE]
                         [--track-high-coverage FILE] [--track-gaps FILE]
                         [--track-min-size TRACK_MIN_SIZE]
                         reference sample

    StrainGR: strain-aware variant caller for metagenomic samples

    positional arguments:
      reference             Reference FASTA file. Can be GZIP compressed.
      sample                BAM file with the aligned reads of the sample against
                            the reference

    optional arguments:
      -h, --help            show this help message and exit

    Quality control:
      Options which determine which reads to consider, when base or read mapping qualities are high enough for calling, etc.

      -Q MIN_QUAL, --min-qual MIN_QUAL
                            Minimum quality for a base to be considered. Default:
                            5.
      -P MIN_PILEUP_QUAL, --min-pileup-qual MIN_PILEUP_QUAL
                            Minimum sum of qualities for an allele to be
                            trusted.for variant calling. Default: 50.
      -F MIN_QUAL_FRAC, --min-qual-frac MIN_QUAL_FRAC
                            Minimum fraction of the reads in the pileup required
                            to confirm an allele (fractions are base quality
                            weighted). Default: 0.1
      -M MIN_MAPPING_QUAL, --min-mapping-qual MIN_MAPPING_QUAL
                            Minimum mapping quality of the whole read to be
                            considered. Default: 5.
      -G MIN_GAP, --min-gap MIN_GAP
                            Minimum size of gap to be considered as such. Default:
                            2000. Will be automatically scaled depending on
                            coverage.

    Output formats:
      Options for writing the results to different file formats. If you're unsure what to choose, output at least the data to HDF5, the rest of the output files can be created afterwards from the HDF5 data using 'straingr view'.

      -s FILE, --summary FILE
                            Output a TSV with a summary of variant calling
                            statistics to the given file.
      --hdf5-out FILE       Output StrainGR variant calling data to the given HDF5
                            file.
      -V FILE, --vcf FILE   Output a VCF file with SNP's. Please be aware that we
                            do not have a good insertion/deletion calling
                            mechanism, but some information on possible indels is
                            written to the VCF file.
      --verbose-vcf         To be used with --vcf. If you set this flag, then the
                            VCF will also include records for positions in the
                            genome where nothing but the reference base is
                            observed. By default it will only output records for
                            positions where some evidence for a SNP is observed.
      --track-callable FILE
                            Output a BED file to the given filename, indicating
                            the callable regions in the genome.
      --track-poor-mq FILE  Output a BED file to the given filename, indicating
                            regions with a majority of poorly mapped reads.
      --track-high-coverage FILE
                            Output a BED file to the given filename, indicating
                            regions with abnormally high coverage.
      --track-gaps FILE     Output a BED file to the given filename, indicating
                            regions marked as a gap
      --track-min-size TRACK_MIN_SIZE
                            For all --track-* options above, only include features
                            (regions) of at least the given size. Default: 1.


#### `strainge sample-compare` - compare variant calls across multiple samples

    usage: strainge sample-compare [-h] [-o SUMMARY_OUT] [-d DETAILS_OUT]
                                   [-b BASELINE] [-D OUTPUT_DIR]
                                   REF SAMPLE_HDF5 [SAMPLE_HDF5 ...]

    Compare strains and variant calls in two different samples. Reads of
    both samples must be aligned to the same reference.

    It's possible to generate a TSV with summary stats as well as a file
    with more detailed information on which alleles are called at what
    positions.

    positional arguments:
      REF                   The reference used for variant calling.
      SAMPLE_HDF5           HDF5 files with variant calling data for each sample.
                            Number of samples should be exactly two, except when
                            used with --baseline.

    optional arguments:
      -h, --help            show this help message and exit
      -o SUMMARY_OUT, --summary-out SUMMARY_OUT
                            Output file for summary statistics. Defaults to
                            stdout.
      -d DETAILS_OUT, --details-out DETAILS_OUT
                            Output file for detailed base level differences
                            between samples (optional).
      -b BASELINE, --baseline BASELINE
                            Path to a sample to use as baseline, and compare all
                            other given samples to this one. Outputs a shell
                            script that runs all individual pairwise comparisons.
      -D OUTPUT_DIR, --output-dir OUTPUT_DIR
                            The output directory of all comparison files when
                            using --baseline.

### Utilities

#### `strainge view` - view variant calling statistics in several file formats

    usage: strainge view [-h] [-s FILE] [--track-covered FILE]
                         [--track-poor-mq FILE] [--track-high-coverage FILE]
                         [--track-gaps FILE] [--track-min-size TRACK_MIN_SIZE]
                         [-V FILE] [--verbose-vcf]
                         reference hdf5

    View call statistics stored in a HDF5 file and output results to
    different file formats

    positional arguments:
      reference             Reference FASTA file. Can be GZIP compressed.
      hdf5                  HDF5 file with StrainGR call statistics.

    optional arguments:
      -h, --help            show this help message and exit
      -s FILE, --summary FILE
                            Output a TSV with a summary of variant calling
                            statistics to the given file.
      --track-covered FILE  Output a BED file to the given filename, indicating
                            callable regions in the genome.
      --track-poor-mq FILE  Output a BED file to the given filename, indicating
                            regions with a majority of poorly mapped reads.
      --track-high-coverage FILE
                            Output a BED file to the given filename, indicating
                            regions marked as abnormally high coverage.
      --track-gaps FILE     Output a BED file to the given filename, indicating
                            regions marked as a gap
      --track-min-size TRACK_MIN_SIZE
                            For all --track-* options above, only include features
                            (regions) of at least the given size. Default: 1.
      -V FILE, --vcf FILE   Output a VCF file with SNP's. Please be aware that we
                            do not have a good insertion/deletion calling
                            mechanism, but some information on possible indels is
                            written to the VCF file.
      --verbose-vcf         To be used with --vcf. If you set this flag, then the
                            VCF will also include rows for positions in the genome
                            where nothing but the reference base is observed. By
                            default it will only output rows for positions where
                            some evidence for a SNP is observed.

#### `strainge stats` - output statistics about k-mer sets

    usage: strainge stats [-h] [-k] [-c] [-H] [-e] [-o OUTPUT] kmerset

    Obtain statistics about a given k-mer set.

    positional arguments:
      kmerset               The K-mer set to load

    optional arguments:
      -h, --help            show this help message and exit
      -k                    Output k-mer size.
      -c, --counts          Output the list of k-mers in this set with
                            corresponding counts.
      -H, --histogram       Write the k-mer frequency histogram to output.
      -e, --entropy         Calculate Shannon entropy in bases and write to
                            output.
      -o OUTPUT, --output OUTPUT
                            Output file, defaults to standard output.

#### `strainge plot` - plot k-mer spectra

    usage: strainge plot [-h] [-o OUTPUT] [-t {spectrum}] kmerset

    Generate plots for a given k-mer set.

    positional arguments:
      kmerset               The k-mer set to load

    optional arguments:
      -h, --help            show this help message and exit
      -o OUTPUT, --output OUTPUT
                            Output filename (PNG preferred).
      -t {spectrum}, --plot-type {spectrum}
                            The kind of plot to generate.
