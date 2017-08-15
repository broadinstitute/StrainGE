#!/usr/bin/env python
"""Generate StrainGR database files

Input: list of fasta genome files
Output: kmer hdf5 files, Bowtie2 index files, and a deduped kmer tree
"""
import os
import sys
import subprocess
import argparse


def run_kmerseq(fasta, k=23, fraction=0.002, force=False):
    """Generate kmer hdf5 file from fasta file"""
    try:
        (root, ext) = os.path.splitext(fasta)
        out = "{}.hdf5".format(root)
        if not force and os.path.isfile(out):
            print >>sys.stderr, "{} already exists, not overwriting".format(out)
            return out
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



def _bowtie2_index_exists(name):
    """Returns true if bowtie2 index exists"""
    for ext in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
        if not os.path.isfile("{}{}".format(name, ext)):
            return False
    return True


def run_bowtie2_build(fasta, force=False):
    """Run bowtie2-build on fasta file"""
    try:
        bowtie2_build = ["bowtie2-build", fasta, fasta]
        if not force and _bowtie2_index_exists(fasta):
            print >>sys.stderr, "Bowtie2 index exists already for {}".format(fasta)
            return True
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


def run_kmertree(kmerfiles, k=23, fingerprint=False, force=False):
    """Generate kmer tree from hdf5 files"""
    try:
        if not force and os.path.isfile("tree.hdf5"):
            print >>sys.stderr, "Warning! Found previously generated kmer tree. Not overwriting"
            return True
        kmertree = ["kmertree", "--dedupe", "-k", str(k), "--output", "tree.hdf5", "--nwk", "tree.nwk"]
        if fingerprint:
            kmertree.append("--fingerprint")
        kmertree.extend(kmerfiles)
        with open("kmertree.log", 'wb') as w:
            subprocess.check_call(kmertree, stdout=w, stderr=w)
        print >>sys.stderr, "Successfully built kmertree"
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
parser.add_argument("-f", "--force", help="Force overwriting database files",
                    action="store_true")
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
        kmerfile = run_kmerseq(fasta, k=args.k, fraction=args.fraction, force=args.force)
        if not kmerfile:
            continue
        kmerfiles.append(kmerfile)

        print >>sys.stderr, "Generating Bowtie2 index for {}".format(fasta)
        if run_bowtie2_build(fasta, force=args.force):
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
    print >>sys.stderr, "ERROR! No kmer files to generate kmer tree"
    sys.exit(1)

if not run_kmertree(kmerfiles, k=args.k, fingerprint=args.fingerprint):
    sys.exit(1)