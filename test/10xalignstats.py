#!/usr/bin/env python
import sys
import pysam

bam = pysam.AlignmentFile(sys.argv[1], "r")
barcodes = {}

def barcode(read):
    for tag, value in read.tags:
        if tag == 'BX':
            return value
    return None


for read in bam:
    if not read.is_unmapped:
        bc = barcode(read)
        bcrange = barcodes.get(bc)
        if bcrange:
            start = min(read.reference_start, bcrange[0])
            end = max(read.reference_end, bcrange[1])
            print bc, barcodes[bc], (start, end)
            barcodes[bc] = (start, end)
        else:
            barcodes[bc] = (read.reference_start, read.reference_end)
            print bc, barcodes[bc]
        #print bc, barcodes[bc]

for bc in barcodes:
    bcrange = barcodes[bc]
    print bc, bcrange[1] - bcrange[0]
