#!/usr/bin/env python

import argparse
import re
import sys
import gzip

import kmertools
import pysam

parser = argparse.ArgumentParser()
parser.add_argument("fastq", help="input fastq with BX barcode tags")
parser.add_argument("input_bam",  help="input bam")
parser.add_argument("output_bam", help="output bam")
args = parser.parse_args()

input = pysam.AlignmentFile(args.input_bam, "r")
output = pysam.AlignmentFile(args.output_bam, "wb", template=input)

tagre = re.compile("BX:Z:([ACGT0-9-]+)")

barcodes = {}

print "Scanning input BAM",
count = 0
with pysam.AlignmentFile(args.input_bam, "r") as bam:
    for read in bam:
        if not read.has_tag("BX"):
            barcodes[read.query_name] = None
        else:
            print read
        count += 1
        if count % 1000000 == 0:
            print "...%d" % (count,),
            sys.stdout.flush()
print
print len(barcodes), "reads need barcodes"

print "Scanning input FASTQ",
count = 0
with gzip.GzipFile(args.fastq, 'r') as fastq:
    for line in fastq:
        if count % 4 == 0:
            assert line[0] == '@', line
            fields = line[1:].strip().split()
            if len(fields) > 1 and fields[1].startswith("BX:Z:") and fields[0] in barcodes:
                barcode = fields[1][5:]
                barcodes[fields[0]] = barcode
        count += 1
        if count % 1000000 == 0:
            print "...%d" % (count,),
            sys.stdout.flush()
fastq.close()
print
print "Found", sum([x is not None for x in barcodes.values()]), '/', len(barcodes), 'barcodes'

print "Writing output BAM",
count = 0
with pysam.AlignmentFile(args.input_bam, "r") as input:
    with pysam.AlignmentFile(args.output_bam, "wb", template=input) as output:
        for read in input:
            name = read.query_name
            if barcodes.get(name, None) and not read.has_tag("BX"):
                read.tags += [("BX", barcodes[name])]
            output.write(read)
            count += 1
            if count % 1000000 == 0:
                print "...%d" % (count,),
                sys.stdout.flush()
print
print 'Done'