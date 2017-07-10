#!/usr/bin/env python


import sys
import pysam

goodbarcodes = []
poorbarcodes = []

# pass 1: collect barcodes with good alignments
with pysam.AlignmentFile(sys.argv[1], "rb") as bam:
    for read in bam:
        barcode = read.get_tag("BX") if read.has_tag("BX") else None
        if barcode and read.has_tag("NM") and read.has_tag("AS"):
            pos = read.get_reference_positions(full_length=True)
            #minpos = min(pos)
            #maxpos = max(pos)
            #conserved = maxpos > 800 and maxpos < 1000 or minpos > 800 and minpos < 1000
            conserved = False
            #mismatches = read.get_tag("NM")
            score = read.get_tag("AS")
            if score > -10 and not conserved:
                goodbarcodes.append(barcode)
            elif score <= -10:
                poorbarcodes.append(barcode)

# pass 2: copy reads with good alignment barcodes
with pysam.AlignmentFile(sys.argv[1], "rb") as bam:
    with pysam.AlignmentFile(sys.argv[2], "wb", template=bam) as out:
        for read in bam:
            barcode = read.get_tag("BX") if read.has_tag("BX") else None
            if barcode in goodbarcodes and barcode not in poorbarcodes:
                out.write(read)


#    output = pysam.AlignmentFile(sys.argv, "wb", template=bt2sam)
