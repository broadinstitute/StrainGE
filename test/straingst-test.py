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
import re

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
    return load_first_column(clusters)


def sample_closest_db(sample_ref, db_refs, similarities):
    if sample_ref in db_refs:
        return sample_ref
    for ref, sim in similarities[sample_ref]:
        if ref in db_refs:
            return ref
    return None


def load_straingst_results(tsvfile):
    with open(tsvfile, 'r') as tsv:
        tsv.__next__()
        tsv.__next__()
        reader = csv.DictReader(tsv, delimiter="\t")
        resultlines = [line for line in reader]
    return resultlines


def sample_index_from_name(filename):
    regex = re.compile(r'sample([0-9]+)')
    return [int(m)-1 for m in regex.findall(filename)]


def eval_stats(truth, results, sample_refs, sample_closest):
    truth_strains = set([sample_refs[t] for t in truth])
    closest_strains = set([sample_closest[t] for t in truth])
    strains = set([r['strain'] for r in results])
    print(truth, truth_strains, closest_strains, strains)
    return len(truth_strains), len(strains)


def analyze(tsvfiles):
    sref = load_sample_refs()
    print(sref)
    dref = load_db_refs()
    sims = load_similarities()
    closest = [sample_closest_db(s, dref, sims) for s in sref]
    print(closest)
    for tsv in glob.glob(tsvfiles):
        truth = sample_index_from_name(tsv)
        results = load_straingst_results(tsv)
        stats = eval_stats(truth, results, sref, closest)
