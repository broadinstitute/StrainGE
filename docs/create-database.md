## StrainGST database creation

### 1. Download high quality reference genomes for your genus/species of interest

This tutorial assumes you have activated the *strainge* conda environment.
The first step is to obtain high quality reference genomes for your genus or
species of interest, any method suffices. We've found the tool [ncbi-genome
-download](https://github.com/kblin/ncbi-genome-download) useful, and will
use that tool for this step. 

For example, to download all *Escherichia* genomes:

```bash
mkdir ref_genomes
ncbi-genome-download bacteria -l complete -g Escherichia,Shigella -H -F all \ 
    -o ref_genomes
```

The `-H` flag automatically organizes all downloaded files in a nice human
-readable folder structure. Besides downloading references, this command
downloads all associated metadata like gene annotations too, which is useful
for downstram analyses.

Next, we organize all references in a single directory using a script
available in the `bin/` directory of this repository: 
`prepare_strainge_db.py`. This script serves two main purposes: 1) it
organizes all references in a single directory, 2) it optionally splits
chromosomes and plasmids into separate files. When tracking strains we're
usually more interested in tracking the chromosome, and we don't want
StrainGST to report a strain as present because it shares a plasmid
(although its algorithm should already prevent most of those cases.)

So download the [prepare_strainge_db.py](https://github.com/broadinstitute/StrainGE/blob/master/bin/prepare_strainge_db.py)
script to your analysis folder, and run it as follows:
```bash
mkdir strainge_db
python3 prepare_strainge_db.py ref_genomes/human_readable -s \
    -o strainge_db > strainge_db/references_meta.tsv
```

The `-s` flag enables splitting chromosomes and plasmids. The file
`references_meta.tsv` contains metadata on each reference (for example its
accession no.)

### 2. K-merize your reference sequences

Next, we k-merize each genome:

```bash
for f in strainge_db/*.fa.gz; do straingst kmerize -o $f.hdf5 $f; done;
```

The FASTA files with only chromosomes have a suffix of `*.chrom.fna.gz`. For
each FASTA file, there is now an accompanying HDF5 file containing the
k-mer data. With `-k` you can optionally specify a different k-mer size, which
by default is 23.

### 3. Compare the k-mer sets and cluster similar references

The goal of StrainGST is to identify close reference genomes to strains present
in a sample. These reference genomes are in turn used for variant calling and
sample comparisons. Here lies a trade-off: the reference genome should be close
enough for accurate variant calling, but sample comparisons are more easy to
perform when the variant calling step is done using the same reference genome,
so you don't want to be too specific. Furthermore, limiting the database
size reduces computational time. The database of reference genomes should 
cover the diversity of the species of interest but not contain too many highly 
similar genomes. Therefore a clustering step is performed to reduce redundancy 
in the database.

We remove redundant reference genomes two ways:

1. Remove reference genomes that are a near perfect subset of another genome.
   An example of this is an *E. coli* strain used for synthetic biology
   applications that was basically a K-12 strain with many genes removed.
2. Cluster closely related genomes based on k-mer similarity and pick one 
   representative.

To do this, we need to compute the pairwise similarities between k-mer sets,
and a metric to identify whether a k-mer set is a subset of another. Both can
be obtained using `straingst kmersim`.

```bash
straingst kmersim --all-vs-all -t 4 -S jaccard -S subset strainge_db/*.hdf5 > similarities.tsv
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
straingst cluster -i similarities.tsv -d -C 0.99 -c 0.90 \
   --clusters-out clusters.tsv \
    strainge_db/*.hdf5 > references_to_keep.txt
```

The cluster command reads our previously created file `similarities.tsv` to
determine which references to keep. The first step is to discard any genome
where more than 99% of its kmers are present in another genome, as enabled by
`-d` and `-C 0.99`. Afterwards, we cluster similar genomes based on the
*Jaccard* similarity between k-mersets: if the Jaccard similarity between two
k-mer sets is higher than 0.90 (`-c 0.90`), those two genomes will be clustered 
together. For each cluster we pick one representative genome: the genome with 
the smallest mean distance to the other cluster members. Each genome to keep is
written to `references_to_keep.txt`. With the option `--clusters-out` we 
specify another file where we write the clustering results. Each line in this 
file specifies a cluster along with its entries, separated by a tab. The
genomes in the first column represent the cluster representatives. This option 
is optional, but can be useful for debugging purposes.

### 4. Create pan-genome k-mer database

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
straingst createdb -o pan-genome-db.hdf5 ref1.hdf5 ref2.hdf5 ...
```

Combining the two methods described above works too.
