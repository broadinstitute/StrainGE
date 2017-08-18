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


def run_kmerseq(fasta, database=None, k=23, fraction=0.002, force=False):
    """Generate kmer hdf5 file from fasta file"""
    if not (fasta and os.path.isfile(fasta)):
        return
    try:
        with open(fasta, 'rb'):
            pass
        if database:
            root = os.path.join(database, os.path.splitext(os.path.basename(fasta)))
            out = os.path.join(database, "{}.hdf5".format(os.path.basename(fasta)))
        else:
            (root, ext) = os.path.splitext(fasta)
            out = "{}.hdf5".format(fasta)
        if not force and os.path.isfile(out):
            print >>sys.stderr, "{} already exists, not overwriting".format(out)
            return out
        kmerseq = ["kmerseq", "-k", str(k), "-o", out, "-f", "--fraction", "{:f}".format(fraction), fasta]
        with open("{}.log".format(root), 'wb') as w:
            subprocess.check_call(kmerseq, stdout=w, stderr=w)
        return out
    except (KeyboardInterrupt, SystemExit) as e:
        raise e
    except IOError:
        print >>sys.stderr, "Cannot open file to kmerize {}".format(fasta)
    except Exception as e:
        print >>sys.stderr, "Exception kmerizing {}:".format(fasta), e


def __kmerseq(cmd):
    """Wrapper to multithread run_kmerseq"""
    try:
        return (cmd[0], run_kmerseq(*cmd))
    except Exception as e:
        raise e



def kmerize_files(fastas, database=None, k=23, fraction=0.002, force=False, threads=1):
    """Kmerize a list of fasta files"""
    try:
        kmerfiles = {}
        if threads > 1:
            p = multiprocessing.Pool(threads)
            cmds = [(fasta, database, k, fraction, force) for fasta in fastas]
            print >>sys.stderr, "Generating kmerized files. Please wait..."
            map_async = p.map_async(__kmerseq, cmds)
            results = map_async.get()
            for (fasta, kmerfile) in results:
                if kmerfile:
                    kmerfiles[kmerfile] = fasta
            p.close()
        else:
            for fasta in fastas:
                try:
                    print >>sys.stderr, "Generating kmerized file for {}".format(fasta)
                    kmerfile = run_kmerseq(fasta, database=database, k=k, fraction=fraction, force=force)
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
        for line in f:
            if line[0] == ">":
                n += 1
    
    return n


def _cluster_kmersim(kmersim, cutoff=0.95):
    """Cluster genomes based on kmer similarity"""
    clusters = {}
    seen = {}
    keep = set()
    root = os.path.split(kmersim)[0]
    with open(kmersim, 'rb') as f:
        print >>sys.stderr, "Clustering kmer similarity results..."
        i = 0
        for line in f:
            temp = line.strip().split("\t")
            if len(temp) != 3:
                continue
            g1 = os.path.join(root, temp[0])
            g2 = os.path.join(root, temp[1])
            keep.update([g1, g2])
            if float(temp[2]) < cutoff:
                continue
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
    
    with open(os.path.join(root, "clusters.txt"), 'wb') as w:
        w.write("\n".join([" ".join(clusters[cluster]) for cluster in clusters]))

    with open(os.path.join(root, "excluded.log"), 'wb') as w:
        for cluster in clusters:
            keep.difference_update(clusters[cluster])
            sort = sorted(clusters[cluster], key = __get_scaffold_count)
            keep.add(sort[0])
            w.write("\n".join(sort[1:])+"\n")
        
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
            print >>sys.stderr, "Running kmer similarity..."
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


def run_bowtie2_build(fasta, database=None, force=False):
    """Run bowtie2-build on fasta file"""
    if not fasta:
        return
    try:
        if database:
            bowtie2_build = ["bowtie2-build", fasta, os.path.join(database, os.path.basename(fasta))]
        else:
            bowtie2_build = ["bowtie2-build", fasta, fasta]
        if not force and _bowtie2_index_exists(fasta):
            print >>sys.stderr, "Bowtie2 index exists already for {}".format(fasta)
            return True
        (root, ext) = os.path.splitext(fasta)
        with open("{}.log".format(root), 'ab') as a:
            subprocess.check_call(bowtie2_build, stdout=a, stderr=a)
        return True
    except (KeyboardInterrupt, SystemExit) as e:
        raise e
    except Exception as e:
        print >>sys.stderr, "Exception building bowtie2 index for {}:".format(fasta), e


def run_kmertree(kmerfiles, k=23, fingerprint=False, force=False):
    """Generate kmer tree from hdf5 files"""
    if not kmerfiles:
        return
    try:
        root = os.path.split(kmerfiles[0])[0]
        treefile = os.path.join(root, "tree.hdf5")
        if not force and os.path.isfile(treefile):
            print >>sys.stderr, "Warning! Found previously generated kmer tree. Not overwriting"
            return True
        kmertree = ["kmertree", "--dedupe", "-k", str(k), "--output", treefile, "--nwk", os.path.join(root, "tree.nwk")]
        if fingerprint:
            kmertree.append("--fingerprint")
        kmertree.extend(kmerfiles)
        with open(os.path.join(root, "kmertree.log"), 'wb') as w:
            print >>sys.stderr, "Generating kmer tree. Please wait..."
            subprocess.check_call(kmertree, stdout=w, stderr=w)
        print >>sys.stderr, "Successfully built kmertree!"
        return True
    except (KeyboardInterrupt, SystemExit):
        print >>sys.stderr, "Interrupting..."
        return
    except Exception as e:
        print >>sys.stderr, "Exception building kmertree:", e


###
### Main
###
def main():
    """Main"""
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", nargs="+", help="reference genome fasta file")
    parser.add_argument("-d", "--database", help="Specify where to save database to (default: in place)")
    parser.add_argument("-k", "--K", help="Kmer size [1-32] (default 23)", type=int, default=23)
    parser.add_argument("--fingerprint", help="use minhash fingerprint instead of full kmer set to build tree (faster for many references)",
                        action="store_true")
    parser.add_argument("--fraction", type=float, default=0.002, help="Fraction of kmers to include in fingerprint (default: 0.002)")
    parser.add_argument("-f", "--force", help="Force overwriting database files",
                        action="store_true")
    parser.add_argument("--cluster", type=float, default=0.95, help="Cluster references at this fraction. Set to 0 to disable (default 0.95)")
    parser.add_argument("--max-contigs", type=int, help="Filter out genomes with more than this many contigs (default: disabled)")
    parser.add_argument("-t", "--threads", type=int, help="Number of threads to use (default: 1)", default=1)
    args = parser.parse_args()

    if not args:
        print >>sys.stderr, "No reference fasta files specified"
        sys.exit(1)

    if args.K > 32:
        print >>sys.stderr, "Cannot use kmer size above 32. Setting to 32..."
        args.K = 32
    elif args.K < 1:
        print >>sys.stderr, "Invalid kmer size of {}. Setting to default of 23...".format(args.K)
        args.K = 23
    
    # record kmer size
    if args.database:
        root = args.database
    else:
        root = os.path.split(args.fasta[0])[0]
    kmersize_file = os.path.join(root, "kmersize")
    if os.path.isfile(kmersize_file):
        with open(kmersize_file, 'rb') as f:
            try:
                kmersize = int(f.readline())
                if kmersize != args.K:
                    print >>sys.stderr, "Overwriting previously generated files with different kmer size of {}".format(kmersize)
                    args.force = True
            except (KeyboardInterrupt, SystemExit):
                print >>sys.stderr, "Interrupting..."
                sys.exit(1)
            except Exception as e:
                print >>sys.stderr, "Exception determining previous kmer size:", e
                print >>sys.stderr, "Will overwrite all previously generated files..."
                args.force = True

    with open(kmersize_file, 'wb') as w:
        w.write("{}\n".format(args.K))

    if args.max_contigs:
        args.fasta = [fasta for fasta in args.fasta if __get_scaffold_count(fasta) <= args.max_contigs]

    complete = 0
    kmerfiles = kmerize_files(args.fasta, database=args.database, k=args.K, fraction=args.fraction, force=args.force, threads=args.threads)
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

    if not run_kmertree(kmerfiles.keys(), k=args.K, fingerprint=args.fingerprint):
        sys.exit(1)

    for ref in kmerfiles:
        try:
            fasta = kmerfiles[ref]
            print >>sys.stderr, "Generating Bowtie2 index for {}".format(fasta)
            if run_bowtie2_build(fasta, database=args.database, force=args.force):
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


if __name__ == "__main__":
    main()