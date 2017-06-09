#!/usr/bin/env python


import sys
import pysam

barcodes = []

# pass 1: collect barcodes with good alignments
with pysam.AlignmentFile(sys.argv[1], "rb") as bam:
    for scaffold in bam:
    for read in bam:
        barcode = read.get_tag("BX") if read.has_tag("BX") else None
        if barcode and read.has_tag("AS") and read.get_tag("AS") > -10:
            barcodes.append(barcode)

# pass 2: copy reads with good alignment barcodes
with pysam.AlignmentFile(sys.argv[1], "rb") as bam:
    with pysam.AlignmentFile(sys.argv[2], "wb", template=bam) as out:
        for read in bam:
            barcode = read.get_tag("BX") if read.has_tag("BX") else None
            if barcode in barcodes:
                out.write(read)


#    output = pysam.AlignmentFile(sys.argv, "wb", template=bt2sam)
