#!/usr/bin/env python

import argparse
import itertools
import re
import subprocess
import sys

import kmertools
import pysam

parser = argparse.ArgumentParser()
parser.add_argument("--fastq1", "-1", required=True, help="read 1 input fastq with BX barcode tags")
parser.add_argument("--fastq2", "-2", required=True, help="read 2 input fastq with BX barcode tags")
parser.add_argument("--reference", "-r", "-x", required=True, help="reference file prefix for bowtie2")
parser.add_argument("--threads", "-p", type=int, default=1, help="number of bowtie2 alignment threads")
parser.add_argument("--output", "-o", help="output bam")
parser.add_argument("--verbose", help="increase output verbosity", action="store_true")
args = parser.parse_args()

BUFSIZE = 16384
tagre = re.compile("BX:Z:([ACGT0-9-]+)")

bowtie2cmd = ["bowtie2", "-p", str(args.threads), "--reorder", "-X", "700", "-x", args.reference,
              "-1", args.fastq1, "-2", args.fastq2, "--reorder"]

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

fastq = kmertools.openSeqFile(args.fastq1)
fq = fastq.next()

for record in bt2sam:
    # ensure we're dealing with the same read name
    if fq.name != record.query_name:
        # step fastq, since it's every other one in sam output
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
    if record.is_paired and record.is_proper_pair:
        aligned_reads += 1
        barcode_aligned = True
    barcode_set.append(record)

# if we have good alignments for this barcode, pass them through
if barcode and barcode_aligned:
    aligned_barcodes += 1
    for r in barcode_set:
        output.write(r)

print aligned_barcodes, "aligned barcodes,", aligned_reads, "aligned reads"
