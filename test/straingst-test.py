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
        resultlines = [line for line in reader if float(line['score']) >= minscore]
    return resultlines


def sample_index_from_name(filename):
    regex = re.compile(r'sample([0-9]+)')
    return [int(m)-1 for m in regex.findall(filename)]


def eval_stats(truth, results, sample_refs, sample_closest):
    truth_strains = set([sample_refs[t] for t in truth])
    closest_strains = set([sample_closest[t] for t in truth])
    strains = set([r['strain'] for r in results])
    db = truth_strains & closest_strains
    tp = strains & closest_strains
    fn = closest_strains - strains
    fp = strains - closest_strains
    return len(tp), len(fn), len(fp), len(db)


def summarize_stats(df, label):
    n = len(df)
    tp = df['TP'].sum()
    fn = df['FN'].sum()
    fp = df['FP'].sum()
    p = tp / (tp + fp)
    r = tp / (tp + fn)
    f = 2 * p * r / (p + r)
    print("{:8s} {:3d} {:3d} {:3d} {:.3f} {:.3f} {:.3f}".format(label, tp, fn, fp, r, p, f))
    return tp, fn, fp, r, p, f


def report_stats(df, label):
    indb_values = df['DB'].unique()
    indb_values.sort()
    for indb in indb_values:
        subset = df[df['DB'] == indb]
        summarize_stats(subset, label + ' ' + str(indb))
    summarize_stats(df, label)


def analyze(tsvfiles, label='All', minscore=0):
    sref = load_sample_refs()
    dref = load_db_refs()
    sims = load_similarities()
    closest = [sample_closest_db(s, dref, sims) for s in sref]
    allstats = []
    for tsv in glob.glob(tsvfiles):
        try:
            truth = sample_index_from_name(tsv)
            results = load_straingst_results(tsv, minscore)
            stats = eval_stats(truth, results, sref, closest)
        except:
            stats = 0, 0, 0, 0
        allstats.append(stats)
    df = pd.DataFrame(allstats, columns=['TP', 'FN', 'FP', 'DB'])
    report_stats(df, label)
    return df

def analyze_all(minscore = 0):
    old_df = None
    total_df = None
    for cov in ['0.1x', '0.5x', '1x', '10x']:
        print("Old")
        old = analyze(os.path.join(cov, '*-' + cov + '-bg.tsv'), label=cov)
        print("New")
        df = analyze(os.path.join(cov, '*-' + cov + '-bg-test.tsv'), minscore=minscore, label=cov)
        old_df = old_df.append(old) if old_df is not None else old
        total_df = total_df.append(df) if total_df is not None else df
    report_stats(old_df, 'Old All')
    report_stats(total_df, 'New All')
    print()


def analyze_minscore(minmin = 0.01, maxmin = 0.02):
    minscore = minmin
    while minscore <= maxmin:
        print("Minscore: " + str(minscore))
        analyze_all(minscore)
        minscore = round(minscore + 0.001, 3)
