#!/usr/bin/env python

import sys
import re
import kmertools
import shelve

read_to_barcode = {}

tagre = re.compile("BX:Z:([ACGT0-9-]+)")

fastq = kmertools.openSeqFile(sys.argv[1])
bcindex = shelve.open("barcode_db.shelf")

count = 0
for fq in fastq:
    id = fq.id
    if id in read_to_barcode:
        print('Dupe read:', id)
        continue
    match = tagre.search(fq.description)
    if match:
        barcode = match.groups()[0]
    else:
        continue
    #print(id, fq.name, barcode)
    bcindex[id] = barcode
    count += 1
    if count % 1000000 == 0:
        print count
bcindex.close()


