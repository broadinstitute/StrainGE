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
    with open(tsvfile, 'r') as tsv:
        tsv.__next__()
        tsv.__next__()
        reader = csv.DictReader(tsv, delimiter="\t")
        resultlines = [line for line in reader if float(line['score']) >= minscore and (len(line['i']) < 3 or line['i'][2] == '0')]
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
            stats = len(tp), len(fn), len(fp), len(db)
        except:
            print("Error reading: " + tsv)
            stats = 0, 0, 0, 0
        if verbose > 2:
            print(f"{stats[0]} {stats[1]} {stats[2]} {stats[3]} {tsv}")
        if fp and verbose > 1:
            tstrains = [sref[t] for t in truth]
            for f in fp:
                dist = [get_ani(t, f, anis) for t in tstrains]
                print(f"{tsv} FP {f} {max(dist)}")
        allstats.append(stats)
    df = pd.DataFrame(allstats, columns=['TP', 'FN', 'FP', 'DB'])
    report_stats(df, label)
    return df


def analyze_1strain(verbose=0, minscore=0, suffix="test", clusters="db/clusters.tsv"):
    old_df = None
    total_df = None
    print("OldNew,Cov,InDB,TP,FN,FP,R,P,F")
    for cov in coverages:
        old = analyze(os.path.join(cov, f"*-{cov}-bg-test"
                                        f".tsv"), label='Old,' + cov, verbose=verbose,
                      clusters="db/clusters.tsv", minscore=minscore)
        df = analyze(os.path.join(cov, f"*-{cov}-bg-{suffix}.tsv"), minscore=minscore,
                     label='New,' + cov, verbose=verbose, clusters=clusters)
        old_df = old_df.append(old) if old_df is not None else old
        total_df = total_df.append(df) if total_df is not None else df
    report_stats(old_df, 'Old,All')
    report_stats(total_df, 'New,All')
    print()


def analyze_minscore(minmin = 0.01, maxmin = 0.02):
    minscore = minminp
    while minscore <= maxmin:
        print("Minscore: " + str(minscore))
        analyze_single(minscore)
        minscore = round(minscore + 0.001, 3)


def analyze_2strain(minscore = 0, verbose = 0, suffix="test", clusters="db/clusters.tsv"):
    total_df = None
    old_df = None
    print("OldNew,Cov,InDB,TP,FN,FP,R,P,F")
    for cov in cov_combinations():
        covdir = f"{cov[0]}-{cov[1]}"
        old = analyze(os.path.join(covdir, "*x-test.tsv"), minscore=minscore, label='Old2S,' + covdir,
                      verbose = verbose, clusters="db/clusters.tsv")
        df = analyze(os.path.join(covdir, f"*-{suffix}.tsv"), minscore=minscore, label='New2S,' + covdir,
                     verbose = verbose, clusters=clusters)
        old_df = old_df.append(old) if old_df is not None else old
        total_df = total_df.append(df) if total_df is not None else df
    report_stats(old_df, 'Old2S,All')
    report_stats(total_df, 'New2S,All')
    print()

def analyze_db(db, minscore = 0):
    dbstr = str(db)
    print("1-strain")
    analyze_1strain(suffix=dbstr, clusters=f"db/clusters{dbstr}.tsv", minscore=minscore)
    #print("2-strain")
    #analyze_2strain(suffix=dbstr, clusters=f"db/clusters{dbstr}.tsv", minscore=minscore)

def db_score_test(db):
    for n in range(10,20):
        s = n / 1000
        print(f"minscore {s}")
        analyze_db(db, minscore=n/1000)
        print()

