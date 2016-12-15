#!/usr/bin/env python

import sys
import urllib

summary = "assembly_summary.txt"

# fetch bacteria assembly_summary.txt file
print("Fetching " + summary)
urllib.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt", summary)

print("Parsing " + summary)
with open(summary, 'r') as s:
    assemblies = [line.strip().split('\t') for line in s]
print(str(len(assemblies)) + " assemblies")
