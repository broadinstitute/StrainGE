#!/usr/bin/env python
import sys
import pysam
import numpy as np

bam = pysam.AlignmentFile(sys.argv[1], "r")
barcodes = {}


histogram = np.zeros((2, 200), dtype=np.int32)

for read in bam:
    if read.is_paired and read.is_proper_pair and read.has_tag("XM"):
        mismatches = read.get_tag("XM")
        which = 0 if read.is_read1 else 1
        histogram[which, mismatches] += 1

for which in xrange(2):
    for count in xrange(200):
        n = histogram[which, count]
        if n:
            print which, count, n
