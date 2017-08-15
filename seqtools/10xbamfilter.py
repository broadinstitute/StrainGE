#!/usr/bin/env python


import sys
import argparse
import pysam


parser = argparse.ArgumentParser()
parser.add_argument("--score", "-s", type=int, default=-10, help="Minimum alignment score")
parser.add_argument("--edge", "-e", type=int, default=250, help="Allow unpaired alignments this close to the edge of contigs")
parser.add_argument("--insert", "-i", type=int, default=-1, help="Minimum insert (fragment) size")
parser.add_argument("--all", "-a", action='store_true', help="Include all matching reads")
parser.add_argument("--verbose", "-v", action='store_true', help="Verbose output")
parser.add_argument("--target", "-t", action='append', help="Targets (e.g., contig3:390-4345)")
parser.add_argument("--exclude", "-x", action='append', help="Targets to exclude(e.g., contig3:390-4345)")
parser.add_argument('input', help='input BAM file')
parser.add_argument('output', help='output BAM file')
args = parser.parse_args()

goodbarcodes = set()
poorbarcodes = set()

def parse_target(target):
    colon = target.rfind(":")
    if colon >= 0:
        contig = target[:colon]
        range_string = target[colon+1:].split("-")
        range_left = int(range_string[0])
        range_right = int(range_string[1]) if len(range_string) > 1 else -1
    else:
        contig = target
        range_left = range_right = -1
    return (contig, range_left, range_right)

def in_target(targets, refname, start, end):
    for t in targets:
        if refname == t[0] and end > t[1] and start < t[2]:
                return t
    return None

includes = map(parse_target, args.target) if args.target else []
excludes = map(parse_target, args.exclude) if args.exclude else []

# pass 1: collect barcodes with good alignments
with pysam.AlignmentFile(args.input, "rb") as bam:
    reflengths = bam.lengths
    refnames = bam.references
    for read in bam:
        barcode = read.get_tag("BX") if read.has_tag("BX") else None
        if barcode and not read.is_unmapped:
            refname = refnames[read.reference_id]
            reflen = reflengths[read.reference_id]
            if includes and not in_target(includes, refname, read.reference_start, read.reference_end):
                continue
            if excludes and in_target(excludes, refname, read.reference_start, read.reference_end):
                continue
            good = True
            score = read.get_tag("AS")
            if score <= args.score:
                good = False
                if args.verbose: print barcode, 'score'
            elif read.is_paired:
                if read.mate_is_unmapped:
                    if read.is_reverse:
                        if read.reference_start > args.edge:
                            good = False
                            if args.verbose: print barcode, 'unmapped mate'
                    else:
                        reflength = reflengths[read.reference_id]
                        if read.reference_end < reflength - args.edge:
                            good = False
                            if args.verbose: print barcode, 'unmapped mate'
                else:
                    if not read.is_proper_pair or abs(read.template_length) < args.insert:
                        good = False
                        if args.verbose: print barcode, 'improper'
            if good:
                goodbarcodes.add(barcode)
                if args.verbose: print barcode, 'good'
            elif not args.all:
                poorbarcodes.add(barcode)

# pass 2: copy reads with good alignment barcodes
with pysam.AlignmentFile(args.input, "rb") as bam:
    with pysam.AlignmentFile(args.output, "wb", template=bam) as out:
        for read in bam:
            barcode = read.get_tag("BX") if read.has_tag("BX") else None
            if barcode in goodbarcodes and barcode not in poorbarcodes:
                out.write(read)


#    output = pysam.AlignmentFile(sys.argv, "wb", template=bt2sam)
