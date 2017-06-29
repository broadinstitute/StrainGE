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
            count, start, end = bcrange
            start = min(read.reference_start, start)
            end = max(read.reference_end, end)
            #print bc, barcodes[bc], (start, end)
            count += 1
            barcodes[bc] = (count, start, end)
        else:
            barcodes[bc] = (1, read.reference_start, read.reference_end)
            #print bc, barcodes[bc]
        #print bc, barcodes[bc]

for bc in barcodes:
    count, start, end = barcodes[bc]
    print bc, count, end - start
