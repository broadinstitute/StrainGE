#  Copyright (c) 2016-2019, Broad Institute, Inc. All rights reserved.
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

import logging
import itertools
from collections import defaultdict

import numpy
import pandas

logger = logging.getLogger(__name__)


def similarities_to_matrix(similarities, labels, metric='jaccard'):
    """Turn the pairwise similarities into a symmetric matrix."""

    label_ix = {label: i for i, label in enumerate(labels)}

    matrix = numpy.empty((len(labels), len(labels)))
    for kmerset1, kmerset2 in similarities.index:
        i = label_ix[kmerset1]
        j = label_ix[kmerset2]

        if metric == 'subset':
            matrix[i, j] = similarities.loc[(kmerset1, kmerset2), 'subset1']
            matrix[j, i] = similarities.loc[(kmerset1, kmerset2), 'subset2']
        else:
            matrix[i, j] = similarities.loc[(kmerset1, kmerset2), metric]
            matrix[j, i] = similarities.loc[(kmerset1, kmerset2), metric]

    for i in range(len(label_ix)):
        matrix[i, i] = numpy.nan

    return pandas.DataFrame(matrix, index=labels, columns=labels)


def cluster_genomes(similarities, labels, threshold, metric='jaccard'):
    label_to_cluster = {
        label: i for i, label in enumerate(labels)
    }

    # Start with a situation where each genome is its own cluster
    clusters = {
        i: [label] for i, label in enumerate(labels)
    }

    # Assumes `similarities` is sorted
    for label1, label2 in similarities.index:
        similarity = similarities.loc[(label1, label2), metric]

        logger.info("Checking %s vs %s, similarity: %g", label1, label2,
                    similarity)

        if similarity < threshold:
            logger.info("Similarity dropped below threshold %g, stopping",
                        threshold)
            break

        if label1 not in label_to_cluster or label2 not in label_to_cluster:
            # Excluded by subset step
            continue

        cluster1 = label_to_cluster[label1]
        cluster2 = label_to_cluster[label2]

        if cluster1 == cluster2:
            # Already in same cluster, skipping
            continue

        clusters[cluster1].extend(clusters[cluster2])

        # Make sure we remember in which cluster each entry lives
        for entry in clusters[cluster2]:
            label_to_cluster[entry] = cluster1

        del clusters[cluster2]

    return clusters


def pick_representative(clusters, similarities, priorities=None,
                        metric='jaccard'):
    if not priorities:
        priorities = {}

    if isinstance(metric, str):
        # Calculate for each label its mean similarity to other cluster members
        sim_per_label = defaultdict(list)
        for cluster, entries in clusters.items():
            if len(entries) == 1:
                sim_per_label[entries[0]].append(1.0)
            else:
                for label1, label2 in itertools.combinations(entries, 2):
                    if (label1, label2) in similarities.index:
                        similarity = similarities.loc[(label1, label2), metric]
                    else:
                        similarity = similarities.loc[(label2, label1), metric]

                    sim_per_label[label1].append(similarity)
                    sim_per_label[label2].append(similarity)

        sort_metric = {
            # Add priority in the tuple to give certain genomes priority when
            # sorting; all else being equal, use shortest label
            label: (priorities.get(label, 1), sum(sim)/len(sim), -len(label))
            for label, sim in sim_per_label.items()
        }
    elif isinstance(metric, dict):
        # It's also possible to specify the metric to sort by directly
        sort_metric = metric
    else:
        raise ValueError(f"Invalid value for parameter 'metric': {metric}")

    # Sort clusters by size (report bigger ones first)
    sorted_clusters = sorted(clusters.items(), key=lambda e: len(e[1]),
                             reverse=True)
    for cluster, entries in sorted_clusters:
        # Report entries with higher mean similarity to other members first
        yield cluster, list(sorted(entries, key=sort_metric.get,
                                   reverse=True))
