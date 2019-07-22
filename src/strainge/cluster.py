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

logger = logging.getLogger(__name__)


def read_similarities(f, threshold=None):
    similarities = {}

    # Assumes similarities file is sorted
    for lineno, line in enumerate(f):
        if line.startswith('#') or not line.strip():
            continue

        parts = [
            p.strip() for p in line.split()
        ]

        if len(parts) != 3:
            raise ValueError(f"Line {lineno} has invalid format. Each line "
                             f"requires three fields.")

        similarity = float(parts[2])

        if threshold and similarity < threshold:
            logger.info("Line %d: Similarity lower than threshold %g",
                        lineno, threshold)
            break

        similarities[tuple(parts[:2])] = float(parts[2])

    return similarities


def cluster_genomes(similarities, labels, threshold):
    label_to_cluster = {
        label: i for i, label in enumerate(labels)
    }

    # Start with a situation where each genome is its own cluster
    clusters = {
        i: [label] for i, label in enumerate(labels)
    }

    # Assumes `similarities` is sorted
    for (label1, label2), similarity in similarities.items():
        logger.info("Checking %s vs %s, similarity: %g", label1, label2,
                    similarity)

        if similarity < threshold:
            logger.info("Similarity dropped below threshold %g, stopping",
                        threshold)
            break

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


def pick_representative(clusters, similarities, priorities=None):
    if not priorities:
        priorities = {}

    # Calculate for each label its mean distance to other cluster members
    sim_per_label = defaultdict(list)
    for cluster, entries in clusters.items():
        if len(entries) == 1:
            sim_per_label[entries[0]].append(1.0)
        else:
            for label1, label2 in itertools.combinations(entries, 2):
                if (label1, label2) in similarities:
                    similarity = similarities[label1, label2]
                else:
                    similarity = similarities[label2, label1]

                sim_per_label[label1].append(similarity)
                sim_per_label[label2].append(similarity)

    entries_mean_dist = {
        # Add priority in the tuple to give certain genomes priority when
        # sorting
        label: (priorities.get(label, 1), sum(sim)/len(sim))
        for label, sim in sim_per_label.items()
    }

    # Sort clusters by size (report bigger ones first)
    sorted_clusters = sorted(clusters.items(), key=lambda e: len(e[1]),
                             reverse=True)
    for cluster, entries in sorted_clusters:
        yield cluster, list(sorted(entries, key=entries_mean_dist.get,
                                   reverse=True))
