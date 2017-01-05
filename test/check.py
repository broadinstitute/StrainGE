#!/usr/bin/env python
import sys
import numpy as np

a = np.load(sys.argv[1])
b = np.load(sys.argv[2])
print a.files, b.files
print 'kmers', all(a['kmers'] == b['kmers'])
print 'counts', all(a['counts'] == b['counts'])

