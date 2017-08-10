#!/usr/bin/env python

import sys
import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

readlist = []

with open(sys.argv[2], "w") as fastq:
    with pysam.AlignmentFile(sys.argv[1], "rb") as bam:
        for read in bam:
            barcode = read.get_tag("BX") if read.has_tag("BX") else ""
            seq = Seq(read.query_sequence)
            if barcode:
                barcode = "BX:Z:" + barcode
            record = SeqRecord(Seq(read.query_sequence), id=read.query_name, description=barcode)
            record.letter_annotations['phred_quality'] = read.query_qualities
            if read.is_reverse:
                record = record.reverse_complement(id=True, description=True)
            readlist.append(record)
    readlist.sort(lambda a, b: cmp(a.description, b.description) or cmp(a.id, b.id) or cmp(len(a), len(b)))
    SeqIO.write(readlist, fastq, "fastq")

for i in xrange(0, len(readlist), 2):
    assert readlist[i].name == readlist[i+1].name
