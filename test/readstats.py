#!/usr/bin/env python


import sys
import pysam
import math
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
import kmertools


parser = argparse.ArgumentParser()
parser.add_argument("--reference", "-r", help="reference file prefix for bowtie2")
parser.add_argument("--bam", "-b", help="bam containing aligned reads")
parser.add_argument("--statsout", "-o", help="hdf5 file in which to save stats")
parser.add_argument("--statsin", "-i", help="hdf5 file from which to load stats")
parser.add_argument("--graph", "-g", help="graph output file (.png best)")
args = parser.parse_args()

READS = 2
POSITIONS = 160
QUALS = 50
MATCHES = 2
# read1/read2, pos, qual, match/mismatch
SHAPE = (READS, POSITIONS, QUALS, MATCHES)
READ_AXIS = 0
POS_AXIS = 1
QUAL_AXIS = 2
MATCH_AXIS = 3


def phred(match, mismatch):
    if match == 0 and mismatch == 0:
        return 0.0
    return round(- 10.0 * math.log10(float(mismatch) / float(mismatch + match)) if mismatch > 0 else 50.0, 1)

def phredprob(q):
    return 10.0 ** (-q/10.0)



def compute_stats(reference_file, bamfile):
    stats = np.zeros(SHAPE, dtype=np.int64)
    print "Loading reference genome"
    # add upper to fix lower case reference genomes
    reference = {scaffold.name: scaffold.seq.upper() for scaffold in kmertools.openSeqFile(reference_file)}
    # pass 1: collect barcodes with good alignments
    with pysam.AlignmentFile(bamfile, "rb") as bam:
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
                if len(column.pileups) > 500:
                    continue
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
                    rpos = pos if not alignment.is_reverse else alignment.query_length - pos - 1
                    qsum += qual
                    if base == refbase:
                        matches.append((read, rpos, qual))
                        matchsum += qual
                    else:
                        mismatches.append((read, rpos, qual))
                if qsum < 100:
                    continue
                matchscore = float(matchsum) / float(qsum)
                if matchscore < .8:
                    continue
                for read, pos, qual in matches:
                    stats[read, pos, qual, 0] += 1
                for read, pos, qual in mismatches:
                    stats[read, pos, qual, 1] += 1
                if refpos % 100000 == 0:
                    print refpos
    return stats

def crunch_stats(stats):
    qsum = stats.sum(axis=QUAL_AXIS)

    fig = plt.figure(1, (5, 10))
    ax1 = fig.add_subplot(211)
    for read in range(READS):
        effqual = []
        calledqual = []
        calledmean = []
        for pos in range(POSITIONS):
            match = 0
            mismatch = 0
            calledprob = 0.0
            calledsum = 0
            for qual in range(QUALS):
                qmatch = stats[read, pos, qual, 0]
                qmismatch = stats[read, pos, qual, 1]
                qtotal = qmatch + qmismatch
                if qtotal:
                    calledprob += phredprob(qual) * qtotal
                    calledsum += qual * qtotal
                    match += qmatch
                    mismatch += qmismatch

            total = match + mismatch
            if not total:
                break
            calledmean.append(float(calledsum) / float(total))
            calledqual.append(phred(total - calledprob, calledprob))
            qual = phred(match, mismatch) if total else 0
            effqual.append(qual)
        ax1.plot(calledmean, label="Read " + str(read + 1) + " mean")
        ax1.plot(calledqual, label="Read " + str(read + 1) + " expected")
        ax1.plot(effqual, label="Read " + str(read + 1) + " actual")
    ax1.legend(loc=3)
    ax1.set_ylim(0, 50)
    ax1.set_ylabel("Base Quality")
    ax1.set_xlabel("Position in Read")

    ax2 = fig.add_subplot(212)
    psum = stats.sum(axis=POS_AXIS)
    for read in range(READS):
        points = []
        for qual in range(QUALS):
            match = float(psum[read, qual, 0])
            mismatch = float(psum[read, qual, 1])
            effqual = phred(match, mismatch) if match + mismatch > 0 else 0
            if effqual:
                points.append((qual, effqual))
        ax2.scatter([p[0] for p in points], [p[1] for p in points], label="Read " + str(read + 1))
    ax2.plot([0, QUALS], [0, QUALS])
    ax2.legend()
    ax1.set_ylim(0, 50)
    ax2.set_ylabel("Effective Base Quality")
    ax2.set_xlabel("Called Base Quality")
    if args.graph:
        plt.savefig(args.graph)
    else:
        plt.show()


############### MAIN

if args.bam and args.reference:
    stats = compute_stats(args.reference, args.bam)
    if args.statsout:
        with h5py.File(args.statsout, 'w') as h5:
            h5.create_dataset("stats", data=stats)
        sys.exit(0)

if args.statsin:
    with h5py.File(args.statsin, 'r') as h5:
        stats = np.array(h5["stats"])
    crunch_stats(stats)





