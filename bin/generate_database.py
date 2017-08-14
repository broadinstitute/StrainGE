#!/usr/bin/env python
"""Generate StrainGR database files

Input: list of fasta genome files
Output: kmer hdf5 files, Bowtie2 index files, and a deduped kmer tree
"""
import os
import sys
import subprocess
import argparse


def run_kmerseq(fasta, k=23, fraction=0.002):
    """Generate kmer hdf5 file from fasta file"""
    try:
        (root, ext) = os.path.splitext(fasta)
        out = "{}.hdf5".format(root)
        kmerseq = ["kmerseq", "-k", str(k), "-o", out, "-f", "--fraction", "{:f}".format(fraction), fasta]
        with open("{}.log".format(root), 'wb') as w:
            subprocess.check_call(kmerseq, stdout=w, stderr=w)
        return out
    except KeyboardInterrupt:
        raise KeyboardInterrupt
    except SystemExit:
        raise SystemExit
    except Exception as e:
        print >>sys.stderr, "Exception kmerizing {}: {}".format(fasta, e)


def run_bowtie2_build(fasta):
    """Run bowtie2-build on fasta file"""
    try:
        bowtie2_build = ["bowtie2-build", fasta, fasta]
        (root, ext) = os.path.splitext(fasta)
        with open("{}.log".format(root), 'ab') as a:
            subprocess.check_call(bowtie2_build, stdout=a, stderr=a)
        return True
    except KeyboardInterrupt:
        raise KeyboardInterrupt
    except SystemExit:
        raise SystemExit
    except Exception as e:
        print >>sys.stderr, "Exception building bowtie2 index for {}: {}".format(fasta, e)


def run_kmertree(kmerfiles, k=23, fingerprint=False):
    """Generate kmer tree from hdf5 files"""
    try:
        if os.path.isfile("tree.hdf5"):
            print >>sys.stderr, "WARNING! Overwriting previously generated kmer tree"
        kmertree = ["kmertree", "--dedupe", "-k", str(k), "--output", "tree.hdf5", "--nwk", "tree.nwk"]
        if fingerprint:
            kmertree.append("--fingerprint")
        kmertree.extend(kmerfiles)
        with open("kmertree.log", 'wb') as w:
            subprocess.check_call(kmertree, stdout=w, stderr=w)
        return True
    except (KeyboardInterrupt, SystemExit):
        print >>sys.stderr, "Interrupting..."
    except Exception as e:
        print >>sys.stderr, "Exception building kmertree: {}".format(e,)


###
### Main
###
parser = argparse.ArgumentParser()
parser.add_argument("fasta", nargs="+", help="reference genome fasta file")
parser.add_argument("-k", "--K", help="Kmer size (default 23)", type=int, default=23)
parser.add_argument("--fingerprint", help="use minhash fingerprint instead of full kmer set to build tree (faster for many references)",
                    action="store_true")
parser.add_argument("--fraction", type=float, default=0.002, help="Fraction of kmers to include in fingerprint (default: 0.002)")
args = parser.parse_args()

if not args:
    print >>sys.stderr, "No reference fasta files specified"
    sys.exit(1)

complete = 0
kmerfiles = []
for fasta in args.fasta:
    try:
        if not os.path.isfile(fasta):
            print >>sys.stderr, "Cannot find file: {}".format(fasta)
            continue
        try:
            open(fasta, 'rb').close()
        except IOError:
            print >>sys.stderr, "Cannot open file: {}".format(fasta)
            continue
        except Exception as e:
            raise e

        print >>sys.stderr, "Generating kmerized file for {}".format(fasta)
        kmerfile = run_kmerseq(fasta, args.database, k=args.k, fraction=args.fraction)
        if not kmerfile:
            continue
        kmerfiles.append(kmerfile)

        print >>sys.stderr, "Generating Bowtie2 index for {}".format(fasta)
        if run_bowtie2_build(fasta):
            complete += 1
    except (KeyboardInterrupt, SystemExit):
        print >>sys.stderr, "Interrupting..."
        sys.exit(1)
    except Exception as e:
        print >>sys.stderr, "Exception processing {}: {}".format(fasta, e)

if complete != len(args.fasta):
    print >>sys.stderr, "ERROR! Processed only {:d} of {:d} references. See log files for details.".format(complete, len(args.fasta))
    sys.exit(1)
else:
    print >>sys.stderr, "Successfully processed all {:d} references".format(complete)

if not kmerfiles:
    print >>sys.stderr, "ERROR! No kmer files to generate kmer tree"`
    sys.exit(1)

if run_kmertree(kmerfiles, k=k, fingerprint=args.fingerprint):
    print >>sys.stderr, "Successfully built kmertree"
else:
    sys.exit(1)