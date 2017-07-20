#!/usr/bin/env python
import sys
import kmertools
import re
from collections import Counter


tag = "BX:Z:"
tagre = re.compile("BX:Z:([ACGT0-9-]+)")


def barcodes(file):
    for seq in kmertools.openSeqFile(file):
        desc = seq.description
        tagindex = desc.find(tag)
        match = tagre.search(desc)
        if match:
            barcode = match.groups()[0]
            yield barcode


counter = Counter(barcodes(sys.argv[1]))
for bc,n in counter.iteritems():
    print bc, n
