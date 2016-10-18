#!/usr/bin/env python
import sys
import kmerizer
import kmertools
import numpy as np
from Bio.Seq import Seq

seq = Seq('agcttttcattctgactgcaacgggcaatatgtctctgtgtggattaaaaaaagagtgtc')
kset = kmertools.KmerSet(int(sys.argv[1]))
kset.kmerizeSeq(str(seq))
print kset.kmers, kset.counts
hashed = kset.hashKmers(kset.kmers)
unhashed = kset.unHashKmers(hashed)
print hashed
print unhashed
print kset.kmers == unhashed
exit(0)
for i in xrange(unhashed.size):
    print "%16x" % kset.kmers[i]
    print "%16x" % hashed[i]
    print "%16x" % unhashed[i]


