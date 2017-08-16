#!/usr/bin/env python
import sys
import argparse

# read wig file for coverage cutoff targets
parser = argparse.ArgumentParser()
parser.add_argument("--cutoff", "-c", type=int, default=100, help="value cutoff")
parser.add_argument("--radius", "-r", type=int, default=200, help="radius around cutoffs")
parser.add_argument("wig", help="Wig file to be parsed")
args = parser.parse_args()
totalcut = 0

def printcut(chrom, start, end):
    global totalcut
    totalcut += end - start + 1
    print "%s:%d-%d" % (chrom, start, end)
    


with open(args.wig, 'r') as wig:
    cutstart = cutend = 0
    coord = 0
    for line in wig:
        tokens = line.strip().split()
        if len(tokens) > 1  and tokens[0] == 'fixedStep':
            if cutstart > 0:
		start = max(cutstart - args.radius, 1)
                printcut(chrom, max(cutstart - args.radius, 1), min(cutend + args.radius, coord - 1))
                cutstart =  cutend = 0
            params = {}
            for t in tokens[1:]:
                thing, value = t.split('=')
                params[thing] = value
            start = int(params.get('start', 1))
            step = int(params.get('step', 1))
            chrom = params.get('chrom')
            coord = start
        elif len(tokens) == 1:
            value = float(tokens[0])
            if value >  args.cutoff:
                cutend = coord
                if cutstart <= 0:
                    cutstart = coord
            else:
                if cutstart > 0 and coord > cutend + args.radius:
                    printcut(chrom, max(cutstart - args.radius, 1), cutend + args.radius)
                    cutstart = cutend = 0
            coord += step


print >>sys.stderr, "Total:", totalcut
