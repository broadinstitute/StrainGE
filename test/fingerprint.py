import sys

import numpy as np
import kmertools

kset = kmertools.KmerSet(int(sys.argv[1]))
kset.load(sys.argv[2])
print kset.fingerprint
for x in kset.fingerprint:
    kset.kmerString(x)


