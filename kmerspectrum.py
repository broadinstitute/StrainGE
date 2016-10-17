#!/usr/bin/env python
import sys
import kmertools
import numpy as np
from optparse import OptionParser

kset = kmertools.KmerSet()
kset.load(sys.argv[1])
print kset.spectrumFilter()
kset.plotSpectrum()

