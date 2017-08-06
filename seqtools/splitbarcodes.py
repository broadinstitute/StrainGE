#!/usr/bin/env python

import os
import sys
import re
import gzip
import shelve
from Bio import SeqIO

def barcode_dir(barcode):
    directory = os.path.join("barcodes", barcode[11:16])
    if not os.path.isdir(directory):
        os.makedirs(directory)
    return directory

def barcode_file(barcode):
    fn = os.path.join(barcode_dir(barcode), barcode + ".fastq.gz")
    if os.path.exists(fn):
        print fn, 'exists!'
    return gzip.GzipFile(fn, "ab")

def write_barcode_files(barcode, reads):
    with barcode_file(barcode) as rfile:
        SeqIO.write(reads, rfile, "fastq")

tagre = re.compile("BX:Z:([ACGT0-9-]+)")
fastq1 = SeqIO.parse(gzip.GzipFile(sys.argv[1], "rb"), "fastq")
if len(sys.argv) > 2:
    fastq2 = SeqIO.parse(gzip.GzipFile(sys.argv[2], "rb"), "fastq")
else:
    fastq2 = fastq1

last_barcode = None
readlist = []

for r1 in fastq1:
    r2 = next(fastq2)
    assert r1.name == r2.name
    match = tagre.search(r1.description)
    if match:
        barcode = match.groups()[0]
        if barcode != last_barcode and readlist:
            print last_barcode
            write_barcode_files(last_barcode, readlist)
            readlist = []
        last_barcode = barcode
        readlist.append(r1)
        readlist.append(r2)
if readlist:
    write_barcode_files(last_barcode, readlist)

    

