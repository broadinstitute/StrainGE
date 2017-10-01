#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


parser = argparse.ArgumentParser()
parser.add_argument("--size", "-s", default=250000, type=int, help="Chunk size")
parser.add_argument("--minsize", "-m", default=5000, type=int, help="Minimum scaffold size")
parser.add_argument("input_fasta", help="input fasta file")
args = parser.parse_args()

for scaffold in SeqIO.parse(open(args.input_fasta, 'r'), "fasta"):
    size = len(scaffold.seq)
    if size < args.minsize: continue
    chunks = (size + args.size - 1) / args.size
    chunksize = (size + chunks - 1) / chunks
    for start in xrange(0, size, chunksize):
        end = min(start + chunksize, size)
        seq = scaffold.seq[start:end]
        name = scaffold.name
        if start > 0 or end < size:
            name += '.' + str(start / chunksize)
        chunk = SeqRecord(seq, id=name, name=name, description="")
        print name, size
        with open(name + ".fasta", 'w') as out:
            SeqIO.write(chunk, out, "fasta")


