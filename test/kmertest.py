#!/usr/bin/env python
import sys
import kmerizer
import kmertools
import numpy as np
from Bio.Seq import Seq

seq = Seq('agcttttcattctgactgcaacgggcaatatgtctctgtgtggattaaaaaaagagtgtc')
k = int(sys.argv[1])
kset = kmertools.KmerSet(k)
kset.kmerizeSeq(str(seq))
print kset.kmers, kset.counts
hashed = kset.hashKmers(kset.kmers)
unhashed = kset.unHashKmers(hashed)
print hashed
print unhashed
print kset.kmers == unhashed

print seq
kmers = kmerizer.kmerize(k, str(seq))
for km in kmers:
    print kmertools.kmerString(k, int(km))
