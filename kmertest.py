#!/usr/bin/env python

import kmerizer
import numpy as np
from Bio.Seq import Seq

seq = Seq('agcttttcattctgactgcaacgggcaatatgtctctgtgtggattaaaaaaagagtgtc')

k15 = kmerizer.kmerize(15, str(seq))

k23 = kmerizer.kmerize(23, str(seq))

k31 = kmerizer.kmerize(31, str(seq))

k23rc = kmerizer.kmerize(23, str(seq.reverse_complement()))

print k23
print k23rc

print k23.size, k23rc.size, np.intersect1d(k23, k23rc).size