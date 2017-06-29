#!/usr/bin/env python

import argparse
import itertools
import re
import subprocess
import sys

import kmertools
import pysam

parser = argparse.ArgumentParser()
parser.add_argument("--fastq1", "-1", help="read 1 input fastq with BX barcode tags")
parser.add_argument("--fastq2", "-2", help="read 2 input fastq with BX barcode tags")
parser.add_argument("--unpaired", "-U", help="read input fastq wiht BX barcode tags")
parser.add_argument("--reference", "-r", "-x", required=True, help="reference file prefix for bowtie2")
parser.add_argument("--threads", "-p", type=int, default=1, help="number of bowtie2 alignment threads")
parser.add_argument("--output", "-o", help="output bam")
parser.add_argument("--verbose", help="increase output verbosity", action="store_true")
args = parser.parse_args()

if args.fastq1 and args.fastq2:
    paired = True
elif args.unpaired:
    paired = False
else:
    parser.print_usage()
    sys.exit(1)

BUFSIZE = 16384
tagre = re.compile("BX:Z:([ACGT0-9-]+)")

if paired:
    bowtie2cmd = ["bowtie2", "-p", str(args.threads), "--reorder", "-X", "700", "-x", args.reference,
                  "-1", args.fastq1, "-2", args.fastq2]
else:
    bowtie2cmd = ["bowtie2", "-p", str(args.threads), "--reorder", "-x", args.reference, "-U", args.unpaired]


bt2 = subprocess.Popen(bowtie2cmd, stdout=subprocess.PIPE, bufsize=BUFSIZE)
bt2sam = pysam.AlignmentFile(bt2.stdout, "r")

if args.output:
    output = pysam.AlignmentFile(args.output, "wb", template=bt2sam)
else:
    output = pysam.AlignmentFile(sys.stdout, "w", template=bt2sam)

barcode = None
barcode_aligned = False
barcode_set = []
aligned_barcodes = 0
aligned_reads = 0

fastq = kmertools.openSeqFile(args.fastq1 if paired else args.unpaired)
fq = fastq.next()

for record in bt2sam:
    # ensure we're dealing with the same read name
    while fq.name != record.query_name:
        fq = fastq.next()
    assert fq.name == record.query_name, "Read names out of sync"

    # extract barcode from input fastq and append it to sam record
    match = tagre.search(fq.description)
    if match:
        new_barcode = match.groups()[0]
        record.tags += [("BX", new_barcode)]
    else:
        new_barcode = None

    if new_barcode != barcode:
        # if we have good alignments for this barcode, pass them through
        if barcode and barcode_aligned:
            aligned_barcodes += 1
            for r in barcode_set:
                output.write(r)
        # start a new set
        barcode = new_barcode
        barcode_set = []
        barcode_aligned = False

    # add this to the barcode set
    if record.mapping_quality > 10 and (record.is_proper_pair or not record.is_paired):
        aligned_reads += 1
        barcode_aligned = True
    barcode_set.append(record)

# if we have good alignments for this barcode, pass them through
if barcode and barcode_aligned:
    aligned_barcodes += 1
    for r in barcode_set:
        output.write(r)

print aligned_barcodes, "aligned barcodes,", aligned_reads, "aligned reads"
