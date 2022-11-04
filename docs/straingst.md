Running StrainGST
-----------------

Identify close reference genome(s) to strain(s) in a sample.

### Prerequisites

1. A pre-built database for the genus or species of interest
2. A whole metagenomic sequencing (WMS) sample

### Usage

#### 1. K-merize the sample reads

StrainGST iteratively compares the k-mer profiles of references in the database
to the k-mers in the sample to identify close reference genomes to strains in a 
sample. 

Our first step is to kmerize the sample reads. For example, if you have a FASTQ 
file named `patient1.fastq`  with all reads, then we generate its corresponding 
k-mer set as follows:

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

We can now run `straingst run` with our database HDF5 and our sample HDF5:

```bash
straingst run -o results.tsv pan-genome-db.hdf5 patient1.hdf5
```

This will output a *tab separated values* (tsv) file, containing statistics 
about the sample k-mer set and a list of identified reference strains with 
accompanying metrics. 

**New in version 1.3:** instead of writing both sample statistics and the identified strains to a
single TSV file, which is generally not as easily read in Python's `pandas` or R, you can now enable 
the option to write sample statistics and strains to separate files when enabling the `--separate-output` (`-O`)
option. If enabling this option, use `-o` to specify the output filename prefix.

Example:

```bash
straingst run -O -o PREFIX pan-genome-df.hdf5 patient1.hdf5
```

This will result in two files: `PREFIX.stats.tsv` (sample statistics), and `PREFIX.strains.tsv` 
(list of identified strains).

### Output file description

#### Example output (single file output)

```
sample    totalkmers      distinct   pkmers  pkcov   pan%
UMB11_01  2277023860      380759656  50090   6.984   1.536
i         strain          gkmers     ikmers  skmers  cov    kcov   gcov   acct   even   spec   rapct  wscore  score
0         Esch_coli_NGF1  49631      49622   50090   0.985  7.009  6.831  0.980  0.987  1.000  1.507  0.940   0.940
```

#### Example output (separate file output; new in version 1.3)

`PREFIX.stats.tsv`

```
sample    totalkmers      distinct   pkmers  pkcov   pan%
UMB11_01  2277023860      380759656  50090   6.984   1.536
```

`PREFIX.strains.tsv`
```
i         strain          gkmers     ikmers  skmers  cov    kcov   gcov   acct   even   spec   rapct  wscore  score
0         Esch_coli_NGF1  49631      49622   50090   0.985  7.009  6.831  0.980  0.987  1.000  1.507  0.940   0.940
```

#### Sample statistics

The first two lines contain statistics on the whole sample. 

Columns:

- *sample*: Sample name, derived from the k-mer set filename
- *totalkmers*: total number of k-mers counted in the sample, including k-mers that occur multiple times
- *distinct*: total _unique_ number of k-mers
- *pkmers*: total unique number of k-mers that are also present in the database
- *pkcov*: average "coverage" (multiplicity) of each unique k-mer in the sample that is also present in the database.
- *pan%*: total number of k-mers (including duplicates) that are present in both the sample and database divided by the
    total number of k-mers in the sample (_totalkmers_), i.e. an estimation of the relative abundance of the
    species/genus of interest in this sample.

#### Reference strain statistics

The next lines contain the close reference genomes identified by StrainGST.

Columns:

- *i*: Iteration number
- *strain*: Reference strain name
- *gkmers*: Total number of unique k-mers in the original reference genome (or its fingerprint). 
- *ikmers*: Remaining unique k-mers in the genome after discarding k-mers excluded in an earlier iteration or because
    their average coverage was too high
- *skmers*: Remaining unique k-mers from the sample
- *cov*: Breadth of coverage of this reference, i.e. what fraction of k-mers in the reference is also present in the
    sample
- *kcov*: Average depth of coverage of k-mers both present in the reference and in the sample
- *gcov*: Average depth of coverage of *all* k-mers in the reference
- *acct*: What fraction of the sample k-mers can be explained by this reference?
- *even*: Evenness of coverage. A value close to 1 indicates that the coverage is evenly distributed along the genome,
    a value close to zero indicates that only a small part of the genome is well covered.
- *spec*: Obsolete
- *rapct*: Estimated strain relative abundance (relative to the whole sample).
- *wscore*: Obsolete
- *score*: Score used to rank each reference in the database at each iteration. A high score represents high confidence
    in this prediction. Scores cannot be compared across iterations or across samples, and it is possible that a strain
    in a second iteration has a higher score than the strain in the first iteration.


### Tips and Tricks

Easily parse StrainGST file in Python (mainly useful for single file output):

```python
from strainge.io.utils import parse_straingst

results = ['sample1.tsv', 'sample2.tsv']

for sample in results:
    print('#', sample)
    with open(sample) as f:
        for strain in parse_straingst(f):
            print(strain)  # strain is a dict with above columns
```

With sample statistics:

```python
from strainge.io.utils import parse_straingst

results = ['sample1.tsv', 'sample2.tsv']

for sample in results:
    print('#', sample)
    with open(sample) as f:
        straingst_iter = iter(parse_straingst(f, return_sample_stats=True))
        sample_stats = next(straingst_iter)
        print(sample_stats)

        for strain in straingst_iter:
            print(strain)  # strain is a dict with above columns
```
