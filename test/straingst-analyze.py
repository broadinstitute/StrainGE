#  Copyright (c) 2016-2020, Broad Institute, Inc. All rights reserved.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
#  * Neither the name Broad Institute, Inc. nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#  POSSIBILITY OF SUCH DAMAGE.
#
import csv
import glob
import os
import re
import pandas as pd
from contextlib import redirect_stdout

#
# Analyzes StrainGST benchmark output
#

coverages = ('10x', '1x', '0.5x', '0.1x')

def cov_combinations():
    return [(coverages[i], coverages[j]) for i in range(len(coverages)) for j in range(i, len(coverages))]


def load_first_column(filename):
    with open(filename, 'r') as f:
        return [line.strip().split()[0] for line in f]


def load_similarities(similarity_file = "db/similarities.tsv"):
    similarities = {}
    with open(similarity_file, 'r') as sim:
        reader = csv.DictReader(sim, delimiter="\t")
        for s in reader:
            s1 = s['kmerset1']
            s2 = s['kmerset2']
            sim = float(s['jaccard'])
            sims1 = similarities.get(s1, [])
            sims1.append((s2, sim))
            similarities[s1] = sims1
            sims2 = similarities.get(s2, [])
            sims2.append((s1, sim))
            similarities[s2] = sims2
    return similarities

def load_anis(similarity_file = "db/similarities.tsv"):
    anis = {}
    with open(similarity_file, 'r') as sim:
        reader = csv.DictReader(sim, delimiter="\t")
        for s in reader:
            s1 = s['kmerset1']
            s2 = s['kmerset2']
            if s1 > s2: s1, s2 = s2, s1
            ani = float(s['ani'])
            anis[(s1, s2)] = ani
    return anis


def get_ani(s1, s2, anis):
    if s1 > s2: s1, s2 = s2, s1
    return anis.get((s1, s2))


def ref_distances():
    refs = list(load_db_refs())
    anis = load_anis()
    results = []
    for i in range(len(refs)):
        for j in range(i + 1, len(refs)):
            s1 = refs[i]
            s2 = refs[j]
            ani = get_ani(s1, s2, anis)
            results.append((ani, s1, s2))
    results.sort()
    return results

def load_sample_refs(ref_list = "samples/ref_list.txt"):
    return load_first_column(ref_list)


def load_db_refs(clusters = "db/clusters.tsv"):
    return set(load_first_column(clusters))


def sample_closest_db(sample_ref, db_refs, similarities):
    if sample_ref in db_refs:
        return sample_ref
    for ref, sim in similarities[sample_ref]:
        if ref in db_refs:
            return ref
    return None


def load_straingst_results(tsvfile, minscore):
    resultlines = []
    with open(tsvfile, 'r') as tsv:
        tsv.__next__()
        tsv.__next__()
        reader = csv.DictReader(tsv, delimiter="\t")
        for line in reader:
            i = line['i']
            if not (len(i) == 1 or i[-1] == '0'):
                continue
            if float(line['score']) < minscore:
                break
            resultlines.append(line)
    return resultlines


def sample_index_from_name(filename):
    regex = re.compile(r'sample([0-9]+)')
    return [int(m)-1 for m in regex.findall(filename)]


def eval_stats(truth, results, sample_refs, sample_closest, verbose=0):
    truth_strains = set([sample_refs[t] for t in truth])
    closest_strains = set([sample_closest[t] for t in truth])
    strains = set([r['strain'] for r in results])
    db = truth_strains & closest_strains
    tp = strains & closest_strains
    fn = closest_strains - strains
    fp = strains - closest_strains
    if verbose > 2 and (fn or fp):
        print(f"Expected: {closest_strains} Found: {strains}")
    return tp, fn, fp, db


def summarize_stats(df, label):
    n = len(df)
    tp = df['TP'].sum()
    fn = df['FN'].sum()
    fp = df['FP'].sum()
    p = tp / (tp + fp)
    r = tp / (tp + fn)
    f = 2 * p * r / (p + r)
    print("{:s},{:d},{:d},{:d},{:f},{:f},{:f}".format(label, tp, fn, fp, r, p, f))
    return tp, fn, fp, r, p, f


def report_stats(df, label):
    indb_values = df['DB'].unique()
    indb_values.sort()
    for indb in indb_values:
        subset = df[df['DB'] == indb]
        summarize_stats(subset, label + ',' + str(indb))
    summarize_stats(df, label + ',All')


def analyze(tsvfiles, label='All', clusters="db/clusters.tsv", minscore=0, verbose=0):
    sref = load_sample_refs()
    dref = load_db_refs(clusters)
    sims = load_similarities()
    anis = load_anis()
    closest = [sample_closest_db(s, dref, sims) for s in sref]
    allstats = []
    for tsv in glob.glob(tsvfiles):
        try:
            truth = sample_index_from_name(tsv)
            results = load_straingst_results(tsv, minscore)
            tp, fn, fp, db = eval_stats(truth, results, sref, closest, verbose=verbose)
            stats = tsv, len(tp), len(fn), len(fp), len(db)
        except:
            print("Error reading: " + tsv)
            stats = tsv, 0, 0, 0, 0
        if verbose > 2:
            print(f"{stats[0]} {stats[1]} {stats[2]} {stats[3]} {tsv}")
        if fp and verbose > 1:
            tstrains = [sref[t] for t in truth]
            for f in fp:
                dist = [get_ani(t, f, anis) for t in tstrains]
                print(f"{tsv} FP {f} {max(dist)}")
        allstats.append(stats)
    df = pd.DataFrame(allstats, columns=['TSV', 'TP', 'FN', 'FP', 'DB'])
    report_stats(df, label)
    return df


def analyze_1strain(suffix="test", verbose=0, minscore=0, clusters="db/clusters.tsv"):
    total_df = None
    print(f"Suffix,Cov,InDB,TP,FN,FP,R,P,F")
    for cov in coverages:
        df = analyze(os.path.join(cov, f"*-{cov}-bg-{suffix}.tsv"), minscore=minscore,
                     label=f"{suffix},{cov}", verbose=verbose, clusters=clusters)
        total_df = total_df.append(df) if total_df is not None else df
    report_stats(total_df, f"{suffix},All")
    return total_df



def analyze_2strain(suffix="test", minscore = 0, verbose = 0, clusters="db/clusters.tsv"):
    total_df = None
    print("Suffix,Cov,InDB,TP,FN,FP,R,P,F")
    for cov in cov_combinations():
        covdir = f"{cov[0]}-{cov[1]}"
        df = analyze(os.path.join(covdir, f"*-{suffix}.tsv"), minscore=minscore, label=f"{suffix},{covdir}",
                     verbose = verbose, clusters=clusters)
        total_df = total_df.append(df) if total_df is not None else df
    report_stats(total_df, f"{suffix},All")
    return total_df



def analyze_db(db, minscore = 0.0, fingerprint=False):
    suffix = dbstr = str(db)
    if fingerprint:
        suffix += "-f"
    if minscore:
        outfile = "analyze-{}-{:.3f}.csv".format(suffix, minscore)
    else:
        outfile = f"analyze-{suffix}.csv"
    with redirect_stdout(open(outfile, 'w')):
        print(f"1-strain minscore={minscore} fingerprint={fingerprint}")
        df1 = analyze_1strain(suffix=suffix, clusters=f"db/clusters{dbstr}.tsv", minscore=minscore)
        print(f"\n2-strain minscore={minscore} fingerprint={fingerprint}")
        df2 = analyze_2strain(suffix=suffix, clusters=f"db/clusters{dbstr}.tsv", minscore=minscore)
    return df1, df2

def db_score_test(db, fingerprint=False):
    for n in range(10,55,5):
        s = n / 1000
        print(f"minscore {s}")
        analyze_db(db, minscore=n/1000, fingerprint=fingerprint)


def dump_closest_refs(db):
    clusters = f"db/clusters{db}.tsv"
    sref = load_sample_refs()
    dref = load_db_refs(clusters)
    sims = load_similarities()
    closest = [sample_closest_db(s, dref, sims) for s in sref]
    with redirect_stdout(open(f"sample-refs-{db}.tsv", 'w')):
        for i, s in enumerate(sref):
            print(f"sample{i+1}\t{s}\t{closest[i]}")
    return sref, closest


def df_to_lookup_dict(df):
    lookup = {}
    pattern = re.compile(".*sample[0-9]+-")
    for i in range(len(df)):
        row = df.iloc[i]
        tsv = row['TSV']
        prefix = pattern.findall(tsv)[0]
        lookup[prefix] = (row['TP'], row['FN'], row['FP'])
    return lookup


def compare_lookup_dict(d1, d2):
    for k in d1:
        if d1[k] != d2[k]:
            print(f"{k}\t{d1[k]}\t{d2[k]}")


def report_misses(db, fingerprint=False, minscore=0.0):
    fp1, fp2 = analyze_db(db, fingerprint=fingerprint, minscore=minscore)
    misses1 = fp1[fp1['FN'] + fp1['FP'] > 0]
    misses2 = fp2[fp2['FN'] + fp2['FP'] > 0]
    misses1.to_csv(f"misses-1strain-{db}-{fingerprint}-{minscore}.csv")
    misses2.to_csv(f"misses-2strain-{db}-{fingerprint}-{minscore}.csv")
    return misses1, misses2


def compare_fingerprint(db, minscore=0):
    full1, full2 = analyze_db(db, fingerprint=False, minscore=minscore)
    full1 = df_to_lookup_dict(full1)
    full2 = df_to_lookup_dict(full2)
    fp1, fp2 = analyze_db(db, fingerprint=True, minscore=minscore)
    fp1 = df_to_lookup_dict(fp1)
    fp2 = df_to_lookup_dict(fp2)
    with redirect_stdout(open(f"fingerprint-diffs-{db}.tsv", 'w')):
        compare_lookup_dict(full1, fp1)
        compare_lookup_dict(full2, fp2)
