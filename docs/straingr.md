Running StrainGR
----------------

Characterize strains in a metagenomic sample.

### Prerequisites

- StrainGST results on one or more samples
- Directory containing the reference genomes used to create the StrainGST database
- BWA-MEM
- Mummer4
- Optional: Picard (for MarkDuplicates)

### Usage

#### 1. Prepare a concatenated reference FASTA with `straingr prepare-ref`

Our strategy to deconvolve strains in a mixture sample is to create a FASTA
containing a close reference genome for each strain present in a sample, and
then aligning the sample reads to this concatenated FASTA file. StrainGR
provides a tool `straingr prepare-ref` to automatically create and analyze
a concatenated reference genome from a list of StrainGST result files.

By including multiple reference genomes into a single FASTA file, reads with an allele 
specific to a strain will be placed to the optimal location in the concatenated reference.
On the other hand, the reference genomes included in the concatenated reference
may share (conserved) parts of their genome because they are the same species, and an 
aligner will be unable to unambiguously place reads in those regions. This is a trade-off: 
include as many reference
genomes as required to deconvolve strains in a sample, without combining too
closely related reference genomes such that they share a vast chunk of their genomes. StrainGR
will not call variants in shared regions.

The `prepare-ref` subcommand aids in building a concatenated reference from
StrainGST result files. It determines which strains have been reported by
StrainGST, and performs another clustering step on the reported strains to
ensure the included reference strains are not too closely related. For 
example, sometimes it happens that a patient has a strain that's somewhat
in the middle between two reference genomes sitting next to each other on the
tree. Due to stochasticity in sequencing, StrainGST may report one reference genome in one
sample, while reporting the other reference in another sample with the same
strain, but taken at a different time point. Here the clustering step ensures
that only one of these two closely related strains gets included in the
concatenated reference.

After concatenating the selected references, `prepare-ref` runs `nucmer` from
the [MUMmer][mummer] toolkit to estimate how "repetitive" the concatenated
reference is, i.e. how much sequence do the genomes concatenated share, by computing
maximal exact matches of at least a configurable size within the concatenated reference.
By default the minimum exact match size is 250 bp, but its recommended to change this value
to the average insert size of the sample read set to most accurately estimate the actual
repetitiveness. These estimates are used to normalize strain abundances in a later step.

[mummer]: https://github.com/mummer4/mummer

To create a concatenated reference, use `straingr prepare-ref` as follows:

```bash
straingr prepare-ref -s path/to/straingst/*.tsv \
   -p "path/to/refdir/{ref}.fa.gz" \
   -S path/to/straingst_db/similarities.tsv
   -o refs_concat.fasta
```

We give multiple StrainGST TSV result files to `prepare-ref` with the `-s`
flag. Usually these are all StrainGST results file belonging to a single
patient, or an other related set of samples. Next, we need to specify how
`prepare-ref` can find the actual FASTA files belonging to strains reported by
StrainGST, this is done using the "path-template" switch `-p`: in this given
path "{ref}" will be replaced with the actual strain name. Don't forgot to use
quotes, because { and } are special characters in many shells. We specify the
similarities.tsv file created at the StrainGST database construction step, to
reuse the calculated k-mer similarities again for clustering. The resulting
concatenated reference will be written to `refs_concat.fasta`.

#### 2. Align reads to the reference

StrainGR is built to be used with `bwa mem`, as it uses the supplied 
information on alternative alignment locations encoded in the `XA` SAM tag to 
deal with shared regions introduced by concatenating reference genomes.

The following command aligns the reads with `bwa mem` and outputs a sorted BAM
file:

```bash
bwa mem -I 300 -t 2 refs_concat.fasta sample1.1.fq.gz sample1.2.fq.gz \
    | samtools sort -@ 2 -O BAM -o sample1.bam -

# Also create BAM index
samtools index sample1.bam
```

We specify a fixed insert size to `bwa mem`, because if the species of interest
in a metagenomic sample is at low abundance, there may be not enough reads per
batch for `bwa mem` to infer the mean insert size, and reads in such a batch 
will be marked as improperly paired. Optionally you can run `picard
MarkDuplicates` on your alignment file.

#### 3. Analyze read alignments to call variants

To call any variants in your sample run the StrainGR variant caller:

```bash
straingr call refs_concat.fasta sample1.bam --hdf5-out sample1.hdf5 --summary sample1.tsv --tracks all
```

All variant calling data will be stored in the given HDF5 file
`sample1.hdf5`. A table with summary statistics like coverage, SNP rate, gaps and more is written to `sample1.tsv`.
You can also specify 
to output this table to a TSV file with the `-s` switch in the above command. 
There are more options for data output, it can output VCF files, BED tracks 
and more, see the CLI reference documentation below.

You can recreate many of the additional data files from the HDF5 file using
`straingr view`.

### Output files

#### StrainGR summary

##### Example output

```
ref                                             name           length    coverage  uReads  abundance  median  callable  callablePct  confirmed  confirmedPct  snps  snpPct  multi  multiPct  lowmq    lowmqPct  high    highPct  gapCount  gapLength
Esch_coli_H3                                    NZ_CP010167.1  4630919   0.247     18      0.000      0       281       0.006        275        97.865        6     2.135   0      0.000     308810   6.668     21045   0.454    13        447553
Esch_coli_H3                                    NZ_CP010168.1  48243     0.025     0       0.000      0       0         0.000        0          0.000         0     0.000   0      0.000     263      0.545     0       0.000    0         0
Esch_coli_NGF1                                  NZ_CP016007.1  5026105   3.549     85824   0.823      3       2506998   49.880       2506921    99.997        77    0.003   70     0.003     1668501  33.197    3681    0.073    1         16868
Esch_coli_NGF1                                  NZ_CP016008.1  40158     6.942     2458    0.008      7       39096     97.355       39094      99.995        2     0.005   2      0.005     982      2.445     12      0.030    0         0
Esch_coli_NGF1                                  NZ_CP016009.1  8556      0.000     0       0.000      0       0         0.000        0          0.000         0     0.000   0      0.000     0        0.000     0       0.000    1         8556
Esch_coli_clone_D_i14                           NC_017652.1    5038386   1.341     210     0.002      0       5022      0.100        5018       99.920        4     0.080   0      0.000     1792289  35.573    196601  3.902    30        575694
Esch_coli_f974b26a-5e81-11e8-bf7f-3c4a9275d6c8  NZ_LR536430.1  4975029   0.224     49      0.000      0       548       0.011        548        100.000       0     0.000   0      0.000     298901   6.008     18373   0.369    16        767577
Esch_coli_1190                                  NZ_CP023386.1  4900891   0.260     24      0.000      0       351       0.007        334        95.157        17    4.843   0      0.000     342936   6.997     24117   0.492    25        848801
Esch_coli_1190                                  NZ_CP023387.1  86147     0.000     0       0.000      0       0         0.000        0          0.000         0     0.000   0      0.000     0        0.000     0       0.000    1         86147
TOTAL                                           -              24754434  1.148     88583   0.834      0       2552296   10.310       2552190    99.996        106   0.004   72     0.003     4412682  17.826    263829  1.066    87        2751196
```

##### Column descriptions

For each scaffold in the concatenated reference StrainGR outputs several metrics.

- *ref*: original reference genome this scaffold belongs to
- *name*: Scaffold name
- *length*: Scaffold length
- *coverage*: Average depth of coverage along this scaffold. **includes multimapped reads, and multimapped reads are
    counted multiple times (for each alternative alignment location)**
- *uReads*: Number of reads uniquely aligned to this scaffold
- *abundance*: Estimated relative abundance of this scaffold. Calculated by dividing the number uniquely aligned reads to
    this scaffold by the total number of reads uniquely aligned, normalized by the estimated repetitiveness in the
    `prepare-ref` stage.
- *median*: median depth of coverage
- *callable* (*callablePct*): Number (percentage) of positions in this scaffold with *strong* evidence for an allele 
    (i.e. two good reads supporting a single allele)
- *confirmed* (*confirmedPct*): Number (percentage) of positions where there's *strong* evidence for the reference 
    allele (does not exclude positions with multiple alleles).
- *snps* (*snpPct*): Number (percentage) of positions with strong evidence for a **single** allele **different** than
    the reference. Our best estimate of ANI.
- *multi* (*multiPct*): Number (percentage) of positions with strong evidence for **multiple alleles** (whether it
    includes the reference or not).
- *lowmq* (*lowmqPct*): Number (percentage) of positions where the majority of reads are mapped with low mapping
    quality, i.e. representing shared or repetitive regions.
- *high* (*highPct*): Number (percentage) of positions with abnormally high coverage.
- *gapCount*: Number of gaps predicted
- *gapLength*: Number of positions in the genome marked as gap


