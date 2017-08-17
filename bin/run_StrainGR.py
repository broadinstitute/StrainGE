#!/usr/bin/env python
"""Run the complete StrainGR Pipeline

Input: fastq file(s)
Output: kmer and genome recovery results
"""
import os
import sys
import subprocess
import multiprocessing

import argparse


def run_kmerseq(fasta, fasta2=None, k=23, fraction=0.002, filtered=False, force=False):
    """Generate kmer hdf5 file from fasta file"""
    try:
        root = os.path.splitext(fasta)[0]
        out = "{}.hdf5".format(fasta)
        if not force and os.path.isfile(out):
            print >>sys.stderr, "{} already exists, not overwriting".format(out)
            return out
        kmerseq = ["kmerseq", "-k", str(k), "-o", out, "-f", "--fraction", "{:f}".format(fraction)]
        if filtered:
            kmerseq.append("-F")
        kmerseq.append(fasta)
        if fasta2:
            kmerseq.append(fasta2)
        with open("{}.log".format(root), 'wb') as w:
            subprocess.check_call(kmerseq, stdout=w, stderr=w)
        return out
    except KeyboardInterrupt:
        raise KeyboardInterrupt
    except SystemExit:
        raise SystemExit
    except Exception as e:
        print >>sys.stderr, "Exception kmerizing {}: {}".format(fasta, e)


def __kmerseq(cmd):
    """Wrapper to multithread run_kmerseq"""
    return (cmd[0], cmd[1], run_kmerseq(*cmd))


def kmerize_files(samples, k=23, fraction=0.002, filtered=False, force=False, threads=1):
    """Kmerize a list of fasta files"""
    try:
        kmerfiles = {}
        if threads > 1 and len(samples) > 1:
            if len(samples) < threads:
                threads = len(samples)
            p = multiprocessing.Pool(threads)
            cmds = []
            for sample in samples:
                if ',' in sample:
                    temp = sample.split(",")
                    if len(temp) != 2:
                        print >>sys.stderr, "ERROR! Could not parse files", sample
                        sys.exit(1)
                    file1, file2 = temp
                else:
                    file1 = sample
                    file2 = None
                cmds.append((file1, file2, k, fraction, filtered, force))
            print >>sys.stderr, "Generating kmerized files. Please wait..."
            map_async = p.map_async(__kmerseq, cmds)
            results = map_async.get()
            for (file1, file2, kmerfile) in results:
                if kmerfile:
                    kmerfiles[kmerfile] = (file1, file2)
            p.close()
        else:
            for sample in samples:
                try:
                    if ',' in sample:
                        temp = sample.split(",")
                        if len(temp) != 2:
                            print >>sys.stderr, "ERROR! Could not parse files", sample
                            sys.exit(1)
                        file1, file2 = temp
                    else:
                        file1 = sample
                        file2 = None

                    print >>sys.stderr, "Kmerizing {}...".format(sample)
                    kmerfile = run_kmerseq(file1, file2, k=k, fraction=fraction, filtered=filtered, force=force)
                    if not kmerfile:
                        continue
                    kmerfiles[kmerfile] = (file1, file2)
                except (KeyboardInterrupt, SystemExit) as e:
                    raise e
                except Exception as e:
                    print >>sys.stderr, "ERROR! Exception while kmerizing {}".format(sample)

        return kmerfiles
    except (KeyboardInterrupt, SystemExit):
        print >>sys.stderr, "Interrupting..."
        if threads > 1 and len(samples) > 1:
            try:
                p.terminate()
            except:
                pass
        sys.exit(1)
    except Exception as e:
        print >>sys.stderr, "Exception while kmerizing files:", e
        if threads > 1 and len(samples) > 1:
            try:
                p.terminate()
            except:
                pass


def run_treepath(kmerfiles, tree, min_score=0.1):
    """Run treepath on a sample kmer file"""
    try:
        treepath = ["treepath", "-o", "treepath.csv", "-s", str(min_score), tree]
        treepath.extend(kmerfiles)
        print >>sys.stderr, "Running path detection. Please wait..."
        with open("treepath.log", 'wb') as w:
            subprocess.check_call(treepath, stdout=w, stderr=w)
        return True
    except (KeyboardInterrupt, SystemExit):
        print >>sys.stderr, "Interrupting..."
    except Exception as e:
        print "ERROR! Exception while running treepath:", e


def parse_treepath():
    """Parse treepath results"""
    results = {}
    if not os.path.isfile("treepath.csv"):
        print >>sys.stderr, "No treepath results found"
        return

    with open("treepath.csv", 'rb') as f:
        f.readline() # skip header
        for line in f:
            temp = line.strip().split(",")
            sample = temp[0]
            strains = []
            for strain in temp[5].split(" "):
                strains.append(":".join(strain.split(":")[:-1]))
            results[sample] = strains
    
    return results


def write_bowtie2_commands(results, kmerfiles, reference, threads=1):
    """Run Bowtie2 aligning samples to references based on treepath"""
    commands = []
    for sample in results:
        file1, file2 = kmerfiles.get(sample)
        for ref in results[sample]:
            index = os.path.join(reference, ref)
            bowtie2 = "bowtie2 --no-unal --very-sensitive --no-mixed --no-discordant -X 700 -p {:d} -x {}".format(threads, index)
            if file2:
                bowtie2 += "-1 {} -2 {}".format(file1, file2)
            else:
                bowtie2 += "-U {}".format(file1)
            bam = "{}_{}.bam".format(sample, ref)
            bowtie2 += " | samtools view -b /dev/stdin | samtools sort -o {};".format(bam)
            bowtie2 += " samtools index {} {}.bai".format(bam, bam)
            commands.append(bowtie2)
    
    with open("bowtie2_commands", 'wb') as w:
        w.write("\n".join(commands))
            

def run_bowtie2(results, kmerfiles, reference, threads=1, force=False):
    """Run Bowtie2 aligner on each sample for each matching reference"""
    total = 0
    aligned = 0
    bamfiles = {}
    for sample in results:
        pair1, pair2 = kmerfiles.get(sample)
        for ref in results[sample]:
            try:
                total += 1
                bam = "{}_{}.bam".format(sample, ref)
                if not force and os.path.isfile(bam):
                    print >>sys.stderr, "Bam file already exists: {}".format(bam)
                    if ref not in bamfiles:
                        bamfiles[ref] = []
                    bamfiles[ref].append(bam)
                    aligned += 1
                    continue
                index = os.path.join(reference, ref)
                bowtie2 = ["bowtie2", "--no-unal", "--very-sensitive", "--no-mixed", "--no-discordant", "-p", str(threads), "-x", index]
                if pair2:
                    print >>sys.stderr, "Aligning {},{} to {}. Please wait...".format(pair1, pair2, ref)
                    bowtie2.extend(["-1", pair1, "-2", pair2])
                else:
                    print >>sys.stderr, "Aligning {} to {}. Please wait...".format(pair1, ref)
                    bowtie2.extend(["-U", pair1])
                
                with open("{}_{}.bowtie2.log".format(sample, ref), 'wb', 0) as w:
                    p_bowtie2 = subprocess.Popen(bowtie2, stdout=subprocess.PIPE, stderr=w)
                    #p_bowtie2.communicate()
                    p_view = subprocess.Popen(["samtools", "view", "-b"], stdin=p_bowtie2.stdout, stdout=subprocess.PIPE, stderr=w)
                    #p_view.communicate()
                    p_sort = subprocess.Popen(["samtools", "sort", "-o", bam], stdin=p_view.stdout, stderr=w)
                    p_sort.communicate()
                    subprocess.check_call(["samtools", "index", bam, "{}.bai".format(bam)])
                    aligned += 1
                    if ref not in bamfiles:
                        bamfiles[ref] = []
                    bamfiles[ref].append(bam)
            except (KeyboardInterrupt, SystemExit):
                print >>sys.stderr, "Interrupting..."
                return
            except Exception as e:
                print "ERROR! Exception occuring during bowtie2 alignment of {} to {}:".format(sample, ref), e
    if aligned == total:
        return bamfiles
    else:
        print >>sys.stderr, "Warning! Only aligned {:d} out of {:d}".format(aligned, total)


def run_straingr(bamfiles, reference):
    """Run straingr tool on bam files to associated reference"""
    complete = 0
    for ref in bamfiles:
        try:
            fasta = os.path.join(reference, ref)
            if not os.path.isfile(fasta):
                print >>sys.stderr, "ERROR! Cannot find reference fasta file: {}".format(fasta)
                continue
            straingr = ["straingr", fasta]
            straingr.extend(bamfiles[ref])
            with open("{}_straingr.log".format(ref), 'wb') as w:
                subprocess.check_call(straingr, stdout=w, stderr=w)
            complete += 1
        except (KeyboardInterrupt, SystemExit):
            print >>sys.stderr, "Interrupting..."
            return
        except Exception as e:
            print >>sys.stderr, "ERROR! Exception occurred during running straingr on {}: {}".format(ref, e)
    
    if complete == len(bamfiles):
        print >>sys.stderr, "Success! All alignments successfully parsed"
        return True
    else:
        print >>sys.stderr, "Error! Not all alignments parsed successfully"





###
### Main
###
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sample", nargs="+", help="sample sequence file(s). If paired end reads, keep pairs together with ','")
    parser.add_argument("-r", "--reference", help="directory containing reference files")
    parser.add_argument("-k", "--K", help="Kmer size (default 23)", type=int, default=23)
    parser.add_argument("-F", "--filter", help="Filter output kmers based on kmer spectrum (to prune sequencing errors at high coverage)",
                        action="store_true")
    parser.add_argument("--fingerprint", help="use minhash fingerprint instead of full kmer set (faster for many references)",
                        action="store_true")
    parser.add_argument("--fraction", type=float, default=0.002, help="Fraction of kmers to include in fingerprint (default: 0.002)")
    parser.add_argument("-s", "--min_score", type=float, default=0.1, help="minimum score of a node in a tree to keep (default: 0.1)")
    parser.add_argument("--no-bowtie2", help="Do not run bowtie2 alignments", 
                        action="store_true")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use")
    parser.add_argument("--force", help="Force overwriting previously generated files",
                        action="store_true")
    args = parser.parse_args()

    if not args.reference:
        print >>sys.stderr, "No reference specified, assuming current working directory"
        args.reference = os.path.curdir

    kmertree = os.path.join(args.reference, "tree.hdf5")
    if not kmertree:
        print >>sys.stderr, "Cannot find tree.hdf5 file in {}".format(args.reference)
        sys.exit(1)

    kmerfiles = kmerize_files(args.sample, k=args.K, fraction=args.fraction, filtered=args.filter, force=args.force, threads=args.threads)
    if not kmerfiles:
        print >>sys.stderr, "ERROR! No kmerized samples found"
        sys.exit(1)

    if not run_treepath(kmerfiles.keys(), kmertree, min_score=args.min_score):
        sys.exit(1)

    treepath_results = parse_treepath()
    if not treepath_results:
        print >>sys.stderr, "ERROR! No treepath results"
        sys.exit(1)

    if args.no_bowtie2:
        write_bowtie2_commands(treepath_results, kmerfiles, args.reference, threads=args.threads)
        print >>sys.stderr, "Wrote bowtie2 alignment commands to bowtie2_commands"
        sys.exit()

    bamfiles = run_bowtie2(treepath_results, kmerfiles, args.reference, threads=args.threads)
    if not bamfiles:
        print >>sys.stderr, "ERROR! Did not complete bowtie2 alignments"
        sys.exit(1)

    if not run_straingr(bamfiles, args.reference):
        sys.exit(1)

    # TODO: downstream analysis!!!

if __name__ == "__main__":
    main()
