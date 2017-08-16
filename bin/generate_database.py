#!/usr/bin/env python
"""Generate StrainGR database files

Input: list of fasta genome files
Output: kmer hdf5 files, Bowtie2 index files, and a deduped kmer tree
"""
import os
import sys
import subprocess
import multiprocessing
import argparse

from Bio import SeqIO


def run_kmerseq(fasta, k=23, fraction=0.002, force=False):
    """Generate kmer hdf5 file from fasta file"""
    if not (fasta and os.path.isfile(fasta)):
        return
    try:
        with open(fasta, 'rb'):
            pass
        (root, ext) = os.path.splitext(fasta)
        out = "{}.hdf5".format(fasta)
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
    except IOError:
        print >>sys.stderr, "Cannot open file to kmerize {}".format(fasta)
    except Exception as e:
        print >>sys.stderr, "Exception kmerizing {}:".format(fasta), e


def __kmerseq(cmd):
    return (cmd[0], run_kmerseq(*cmd))


def kmerize_files(fastas, k=23, fraction=0.002, force=False, threads=1):
    """Kmerize a list of fasta files"""
    try:
        kmerfiles = {}
        if threads > 1:
            p = multiprocessing.Pool(threads)
            cmds = [(fasta, k, fraction, force) for fasta in fastas]
            print >>sys.stderr, "Generating kmerized files. Please wait..."
            map_async = p.map_async(__kmerseq, cmds)
            results = map_async.get()
            for (fasta, kmerfile) in results:
                kmerfiles[kmerfile] = fasta
            p.close()
        else:
            for fasta in args.fasta:
                try:
                    print >>sys.stderr, "Generating kmerized file for {}".format(fasta)
                    kmerfile = run_kmerseq(fasta, k=args.K, fraction=args.fraction, force=args.force)
                    if not kmerfile:
                        continue
                    kmerfiles[kmerfile] = fasta
                except (KeyboardInterrupt, SystemExit) as e:
                    raise e
                except Exception as e:
                    print >>sys.stderr, "Exception kmerizing {}:".format(fasta), e

        return kmerfiles
    except (KeyboardInterrupt, SystemExit):
        print >>sys.stderr, "Interrupting..."
        if threads > 1 and p:
            try:
                p.terminate()
            except:
                pass
        sys.exit(1)
    except Exception as e:
        print >>sys.stderr, "Exception while kmerizing files:", e
        if threads > 1 and p:
            try:
                p.terminate()
            except:
                pass





def __get_scaffold_count(fasta):
    """Get the number of scaffolds in a fasta file"""
    if not fasta:
        print >>sys.stderr, "Cannot find fasta file {}".format(fasta)
        raise IOError
    n = 0
    with open(fasta, 'rU') as f:
        for seq in SeqIO.parse(f, "fasta"):
            n += 1
    
    return n


def _cluster_kmersim(kmersim, cutoff=0.95):
    """Cluster genomes based on kmer similarity"""
    clusters = {}
    seen = {}
    keep = set()
    with open(kmersim, 'rb') as f:
        print >>sys.stderr, "Clustering kmer similarity results"
        i = 0
        for line in f:
            temp = line.strip().split("\t")
            if len(temp) != 3:
                continue
            keep.update(temp[:2])
            if float(temp[2]) < cutoff:
                break # file is sorted...
            g1, g2 = temp[:2]
            keep.remove(g1)
            keep.remove(g2)
            if g1 in seen and g2 not in seen:
                c1 = seen[g1]
                clusters[c1].add(g2)
                seen[g2] = c1
            elif g2 in seen and g1 not in seen:
                c2 = seen[g2]
                clusters[c2].add(g1)
                seen[g1] = c2
            elif g1 in seen and g2 in seen:
                c1 = seen[g1]
                c2 = seen[g2]
                if c1 < c2:
                    clusters[c1].update(clusters[c2])
                    for g in clusters[c2]:
                        seen[g] = c1
                    del clusters[c2]
                elif c2 < c1:
                    clusters[c2].update(clusters[c1])
                    for g in clusters[c1]:
                        seen[g] = c2
                    del clusters[c1]
            else:
                clusters[i] = set([g1, g2])
                seen[g1] = i
                seen[g2] = i
                i += 1
    
    for cluster in clusters:
        keep.add(min(clusters[cluster], key = lambda x: __get_scaffold_count(os.path.splitext(x)[0])))
    
    print >>sys.stderr, "After clustering, {:d} genomes remain".format(len(keep))
    return ["{}.hdf5".format(name) for name in keep]
    


def run_kmersim(kmerfiles, fingerprint=False, threads=1, cutoff=0.95, force=False):
    """Run kmer similarity tool on reference hdf5 files"""
    if not kmerfiles:
        return
    if cutoff <= 0:
        return kmerfiles
    try:
        root = os.path.split(kmerfiles[0])[0]
        out = os.path.join(root, "kmersim.txt")
        if not force and os.path.isfile(out):
            print >>sys.stderr, "kmersimilarity results already exist"
            return _cluster_kmersim(out)
        kmersim = ["kmersimilarity", "--similarity", out]
        if threads > 1:
            kmersim.extend(["--threads", str(threads)])
        if fingerprint:
            kmersim.append("--fingerprint")
        kmersim.extend(kmerfiles)
        with open(os.path.join(root, "kmersim.log"), 'wb') as w:
            subprocess.check_call(kmersim, stdout=w, stderr=w)
        return _cluster_kmersim(out)
    
    except (KeyboardInterrupt, SystemExit):
        print >>sys.stderr, "Interrupting..."
        sys.exit(1)
    except Exception as e:
        print >>sys.stderr, "Exception running kmersimilarity: ", e
    


def _bowtie2_index_exists(name):
    """Returns true if bowtie2 index exists"""
    if not name:
        return
    for ext in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
        if not os.path.isfile("{}{}".format(name, ext)):
            return
    return True


def run_bowtie2_build(fasta, force=False):
    """Run bowtie2-build on fasta file"""
    if not fasta:
        return
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
        print >>sys.stderr, "Exception building bowtie2 index for {}:".format(fasta), e


def run_kmertree(kmerfiles, k=23, fingerprint=False, force=False):
    """Generate kmer tree from hdf5 files"""
    if not kmerfiles:
        return
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
        print >>sys.stderr, "Exception building kmertree:", e


###
### Main
###
def main():
    """Main"""
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", nargs="+", help="reference genome fasta file")
    parser.add_argument("-k", "--K", help="Kmer size (default 23)", type=int, default=23)
    parser.add_argument("--fingerprint", help="use minhash fingerprint instead of full kmer set to build tree (faster for many references)",
                        action="store_true")
    parser.add_argument("--fraction", type=float, default=0.002, help="Fraction of kmers to include in fingerprint (default: 0.002)")
    parser.add_argument("-f", "--force", help="Force overwriting database files",
                        action="store_true")
    parser.add_argument("--cluster", type=float, default=0.95, help="Cluster references at this fraction. Set to 0 to disable (default 0.95)")
    parser.add_argument("-t", "--threads", type=int, help="Number of threads to use (default: 1)", default=1)
    args = parser.parse_args()

    if not args:
        print >>sys.stderr, "No reference fasta files specified"
        sys.exit(1)

    complete = 0
    kmerfiles = kmerize_files(args.fasta, k=args.K, fraction=args.fraction, force=args.force, threads=args.threads)
    if not kmerfiles:
        sys.exit(1)
    if len(kmerfiles) != len(args.fasta):
        print >>sys.stderr, "ERROR! Only kmerized {:d} out of {:d} files".format(len(kmerfiles), len(args.fasta))
        sys.exit(1)

    if args.cluster > 0:
        keep = run_kmersim(kmerfiles.keys(), fingerprint=args.fingerprint, threads=args.threads, cutoff=args.cluster, force=args.force)
        if not keep:
            sys.exit(1)
        for ref in kmerfiles.keys():
            if ref not in keep:
                # remove unselected reference genomes
                del kmerfiles[ref]

    for ref in kmerfiles:
        try:
            fasta = kmerfiles[ref]
            print >>sys.stderr, "Generating Bowtie2 index for {}".format(fasta)
            if run_bowtie2_build(fasta, force=args.force):
                complete += 1

        except (KeyboardInterrupt, SystemExit):
            print >>sys.stderr, "Interrupting..."
            sys.exit(1)
        except Exception as e:
            print >>sys.stderr, "Exception building bowtie2 index for {}: {}".format(fasta, e)

    if complete != len(kmerfiles):
        print >>sys.stderr, "ERROR! Processed only {:d} of {:d} references. See log files for details.".format(complete, len(args.fasta))
        sys.exit(1)
    else:
        print >>sys.stderr, "Successfully processed all {:d} references".format(complete)

    if not kmerfiles:
        print >>sys.stderr, "ERROR! No kmer files to generate kmer tree"
        sys.exit(1)

    if not run_kmertree(kmerfiles.keys(), k=args.k, fingerprint=args.fingerprint):
        sys.exit(1)


if __name__ == "__main__":
    main()