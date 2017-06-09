#!/usr/bin/env python


import sys
import pysam
import math
import kmertools
import numpy as np

print "Loading reference genome"
# add upper to fix lower case reference genomes
reference = {scaffold.name: scaffold.seq.upper() for scaffold in kmertools.openSeqFile(sys.argv[1])}
print reference.keys()

# match/mismatch, read1/read2, pos, qual
stats = np.zeros((2, 160, 50, 2), dtype=np.int64)

# pass 1: collect barcodes with good alignments
with pysam.AlignmentFile(sys.argv[2], "rb") as bam:
    for scaffold in bam.references:
        refseq = reference[scaffold]
        if len(refseq) < 1000000: continue
        print scaffold, len(refseq)
        for column in bam.pileup(scaffold):
            refpos = column.reference_pos
            refbase = refseq[refpos]
            if refbase == 'N':
                continue
            matches = []
            mismatches = []
            qsum = 0
            matchsum = 0
            for read in column.pileups:
                alignment = read.alignment
                if alignment.is_unmapped:
                    continue
                if alignment.is_paired and not alignment.is_proper_pair:
                    continue
                pos = read.query_position_or_next
                qual = alignment.query_qualities[pos]
                base = alignment.query_sequence[pos]
                read = 0 if alignment.is_read1 else 1
                qsum += qual
                if base == refbase:
                    matches.append((read, pos, qual))
                    matchsum += qual
                else:
                    mismatches.append((read, pos, qual))
            if qsum < 100:
                continue
            matchscore = float(matchsum) / float(qsum)
            if matchscore < .75:
                continue
            for read, pos, qual in matches:
                stats[read, pos, qual, 0] += 1
            for read, pos, qual in mismatches:
                stats[read, pos, qual, 1] += 1
            if refpos % 1000 == 0:
                print refpos

for read in range(2):
    for pos in range(160):
        for qual in range(50):
            match = float(stats[read, pos, qual, 0])
            mismatch = float(stats[read, pos, qual, 1])
            total = match + mismatch
            if total > 0:
                effqual = - 10.0 * math.log10(mismatch / total) if mismatch > 0 else 50.0
                effqual = round(effqual, 1)
                print read, pos, qual, effqual




