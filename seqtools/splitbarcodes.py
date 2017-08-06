#!/usr/bin/env python

import sys
import kmertools
import gzip
from Bio import SeqIO

import sys
import kmertools

with gzip.GzipFile(sys.argv[1], "rb") as read1, gzip.GzipFile(sys.argv[2], "rb") as read2:

    lastread = None
    for read in read1:
        if lastread:
            if lastread.name == read.name:
                SeqIO.write(lastread, read1, "fastq")
                SeqIO.write(read, read2, "fastq")
                lastread = None
            else:
                print 'No Match:', lastread
        else:
            lastread = read