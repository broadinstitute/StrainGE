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
parser.add_option("-b", "--batch", type='int', default=100000000, help='Process input kmers in batches of this size')
parser.add_option("-s", "--spectrum", type='str', help='File to store kmer spectrum png')
parser.add_option("-m", "--smax", type='int', default=300, help="Maximum frequency in spectrum plot")
parser.add_option("-o", "--output", type='str', help="Output file for kmers and counts")
parser.add_option("-c", "--compress", action='store_true', help="Store kmer output file in compressed format")

(options, args) = parser.parse_args()


def openFastq(fileName):
    if fileName.endswith(".bz2"):
        return bz2.BZ2File(fileName, 'r')
    elif fileName.endswith(".gz"):
        return gzip.GzipFile(fileName, 'r')
    return open(fileName, 'r')

def mergeCounts(unique, counts, newKmers):
    newUnique, newCounts = np.unique(kmerBatch[:nkmers], return_counts=True)
    print 'New:', newUnique.size, newCounts.sum()
    if type(unique) == type(None):
        return (newUnique, newCounts)
    else:
        maxSize = unique.size + newUnique.size
        mergedUnique = np.zeros(maxSize, dtype=np.uint64)
        mergedCounts = np.zeros(maxSize, dtype=np.int64)
        count = kmerizer.merge_counts(unique, counts, newUnique, newCounts, mergedUnique, mergedCounts)
        print 'Merged:', count, mergedCounts.sum()
        return (mergedUnique[:count], mergedCounts[:count])

kmers = None
counts = None

kmerBatch = np.zeros(options.batch, dtype=np.uint64)

nreads = 0
nbases = 0
nkmers = 0
totalReads = 0
totalBases = 0
totalKmers = 0

def processBatch():
    global kmerBatch, nkmers, nbases, nreads, totalKmers, totalBases, totalReads, kmers, counts
    totalReads += nreads
    totalBases += nbases
    totalKmers += nkmers
    print 'Reads:', totalReads, ', Bases:', totalBases, ', Kmers:', totalKmers
    kmers, counts = mergeCounts(kmers, counts, kmerBatch[:nkmers])
    nreads = 0
    nbases = 0
    nkmers = 0


for arg in args:
    fastq = SeqIO.parse(openFastq(arg), "fastq")
    for read in fastq:
        nkmers += kmerizer.kmerize_into_array(options.k, str(read.seq), kmerBatch, nkmers)
        nreads += 1
        nbases += len(read)
        if nbases > options.batch:
            processBatch()
    fastq.close()
if nkmers:
    processBatch()
kmerBatch = kmers
print 'Unique Kmers:', kmerBatch.size, ', Singletons:', kmerBatch.size - np.count_nonzero(counts - 1)

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

if options.output:
    print 'Writing output to', options.output
    if options.compress:
        np.savez_compressed(options.output, kmers=kmers, counts=counts)
    else:
        np.savez(options.output, kmers=kmers, counts=counts)
