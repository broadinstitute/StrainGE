#!/usr/bin/env python
import sys
import gzip
import bz2
from itertools import combinations
from Bio import SeqIO
import kmerizer
import numpy as np
import pydot
from optparse import OptionParser
import matplotlib.pyplot as plt

parser = OptionParser()
parser.add_option("-k", type='int', help="Kmer size", default=31)
parser.add_option("-s", "--spectrum", type='str', help='File to store kmer spectrum png')
parser.add_option("--smax", type='int', help="Maximum frequency in spectrum plot")

(options, args) = parser.parse_args()


def openFastq(fileName):
    if fileName.endswith(".bz2"):
        return bz2.BZ2File(fileName, 'r')
    elif fileName.endswith(".gz"):
        return gzip.GzipFile(fileName, 'r')
    return open(fileName, 'r')

kmers = np.zeros(1000000000, dtype=np.uint64)
nreads = 0
nbases = 0
nkmers = 0

for arg in args:
    fastq = SeqIO.parse(openFastq(arg), "fastq")
    for read in fastq:
        nkmers += kmerizer.kmerize_into_array(options.k, str(read.seq), kmers, nkmers)
        nreads += 1
        nbases += len(read)
        if nbases > 1000000000:
            break
    fastq.close()
    print 'Reads:', nreads, ', Bases:', nbases, ', Kmers:', nkmers
(kmers, counts) = np.unique(kmers[:nkmers], return_counts=True)
print 'Unique Kmers:', kmers.size, ', Singletons:', kmers.size - np.count_nonzero(counts - 1)

if options.spectrum:
    # to get kmer profile, count the counts!
    profile = np.unique(counts, return_counts=True)
    plt.semilogy(profile[0], profile[1])
    plt.grid = True
    if options.smax:
        plt.xlim(1, options.smax)
    plt.xlabel("Kmer Frequency")
    plt.ylabel("Number of Kmers")
    plt.suptitle = "Kmer Spectrum (K=%d)" % (options.k,)
    plt.savefig(options.spectrum)



