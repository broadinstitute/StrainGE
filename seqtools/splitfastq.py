#!/usr/bin/env python

import sys
import kmertools
import gzip
from Bio import SeqIO

import sys
import kmertools

with open(sys.argv[2], "w") as read1, open(sys.argv[3], "w") as read2:
    lastread = None
    for read in kmertools.openSeqFile(sys.argv[1]):
        if lastread:
            if lastread.name == read.name:
                SeqIO.write(lastread, read1, "fastq")
                SeqIO.write(read, read2, "fastq")
                lastread = None
            else:
                print 'No Match:', lastread
        else:
            lastread = read

