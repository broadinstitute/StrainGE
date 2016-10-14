#!/usr/bin/env python

import kmerizer
import numpy as np
from Bio.Seq import Seq

seq = Seq('agcttttcattctgactgcaacgggcaatatgtctctgtgtggattaaaaaaagagtgtc')
seq2 = seq[25:] + seq[:25]

k15 = kmerizer.kmerize(15, str(seq))

k23 = kmerizer.kmerize(23, str(seq))
k232 = kmerizer.kmerize(23, str(seq2))

k31 = kmerizer.kmerize(31, str(seq))

k23rc = kmerizer.kmerize(23, str(seq.reverse_complement()))

print k23
print k232

print k23.size, k23rc.size, np.intersect1d(k23, k232).size


k23counts = np.unique(k23, return_counts=True)
k232counts = np.unique(k232, return_counts=True)
countSum = k23.size + k232.size
mergedKmers = np.zeros(countSum, dtype=np.uint64)
mergedCounts = np.zeros(countSum, dtype=np.int64)
print k23counts
print k232counts
print mergedKmers, mergedCounts
count = kmerizer.merge_counts(k23counts[0], k23counts[1], k232counts[0], k232counts[1], mergedKmers, mergedCounts)
mergedKmers = mergedKmers[:count]
mergedCounts = mergedCounts[:count]
print count, mergedKmers, mergedCounts, mergedCounts.sum()
