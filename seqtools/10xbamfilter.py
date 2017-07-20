#!/usr/bin/env python


import sys
import argparse
import pysam


parser = argparse.ArgumentParser()
parser.add_argument("--score", "-s", type=int, default=-10, help="Minimum alignment score")
parser.add_argument("--edge", "-e", type=int, default=250, help="Allow unpaired alignments this close to the edge of contigs")
parser.add_argument("--insert", "-i", type=int, default=-1, help="Minimum insert (fragment) size")
parser.add_argument("--contig", "-c", type=int, default=-1, help="Only filter based on this contig (index)")
parser.add_argument('input', help='input BAM file')
parser.add_argument('output', help='output BAM file')
args = parser.parse_args()

goodbarcodes = set()
poorbarcodes = set()

# pass 1: collect barcodes with good alignments
with pysam.AlignmentFile(args.input, "rb") as bam:
    reflengths = bam.lengths
    for read in bam:
        barcode = read.get_tag("BX") if read.has_tag("BX") else None
        if barcode and not read.is_unmapped:
            if args.contig >= 0 and read.reference_id != args.contig:
                continue
            good = True
            score = read.get_tag("AS")
            if score <= args.score:
                good = False
            elif read.is_paired:
                if read.mate_is_unmapped:
                    if read.is_reverse:
                        if read.reference_start > args.edge:
                            good = False
                    else:
                        reflength = reflengths[read.reference_id]
                        if read.reference_end < reflength - args.edge:
                            good = False
                else:
                    if not read.is_proper_pair or abs(read.template_length) < args.insert:
                        good = False
            if good:
                goodbarcodes.add(barcode)
            else:
                poorbarcodes.add(barcode)

# pass 2: copy reads with good alignment barcodes
with pysam.AlignmentFile(args.input, "rb") as bam:
    with pysam.AlignmentFile(args.output, "wb", template=bam) as out:
        for read in bam:
            barcode = read.get_tag("BX") if read.has_tag("BX") else None
            if barcode in goodbarcodes and barcode not in poorbarcodes:
                out.write(read)


#    output = pysam.AlignmentFile(sys.argv, "wb", template=bt2sam)
