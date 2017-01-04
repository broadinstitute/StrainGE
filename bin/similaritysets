#!/usr/bin/env python

import sys
import csv

reader = csv.reader(open(sys.argv[1], 'r'), delimiter="\t")
limit = float(sys.argv[2])

setIndex = {}
sims = []
for line in reader:
    (strain1, strain2, sim) = line
    sim = float(sim)
    if strain1 not in setIndex:
        setIndex[strain1] = {strain1}
    if strain2 not in setIndex:
        setIndex[strain2] = {strain2}
    if sim < limit:
        continue
    sims.append((strain1, strain2, sim))
    # merge sets
    set1 = setIndex[strain1]
    set2 = setIndex[strain2]
    if set1 != set2:
        set1 |= set2
        for s in set2:
            setIndex[s] = set1

def simScore(genome):
    gsims = [sim[2] for sim in sims if sim[0] == genome or sim[1] == genome]
    return sum(gsims)/float(len(gsims))
    
sets = []

for s in setIndex.values():
    if s not in sets:
        sets.append(s)
        if len(s) == 1:
            print " ".join(s)
        else:
            slist = [(g, simScore(g)) for g in s]
            slist.sort(lambda a, b: cmp(b[1], a[1]))
            print " ".join([g[0] for g in slist])
print >>sys.stderr, len(sets), 'sets containing', sum([len(s) for s in sets]), 'genomes'
