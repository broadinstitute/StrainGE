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

import os
import sys
import csv
import logging
import argparse
import functools
import itertools
import multiprocessing

import h5py

from strainge.cli.registry import Subcommand
from strainge import kmertools, utils, comparison

logger = logging.getLogger(__name__)


class StatsSubcommand(Subcommand):
    """
    Obtain statistics about a given k-mer set.
    """

    def register_arguments(self, subparser: argparse.ArgumentParser):
        subparser.add_argument(
            'kmerset',
            help="The K-mer set to load"
        )
        subparser.add_argument(
            '-k', action="store_true", default=False,
            help="Output k-mer size."
        )
        subparser.add_argument(
            '-c', '--counts', action="store_true", default=False,
            help="Output the list of k-mers in this set with corresponding "
                 "counts."
        )
        subparser.add_argument(
            '-H', '--histogram', action="store_true", default=False,
            help="Write the k-mer frequency histogram to output."
        )
        subparser.add_argument(
            '-e', '--entropy', action="store_true", default=False,
            help="Calculate Shannon entropy in bases and write to output."
        )
        subparser.add_argument(
            '-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
            help="Output file, defaults to standard output."
        )

    def __call__(self, kmerset, output, k=False, counts=False, histogram=False,
                 entropy=False, **kwargs):
        logger.info("Loading k-merset %s", kmerset)
        kmerset = kmertools.kmerset_from_hdf5(kmerset)

        if k:
            print("K", kmerset.k, file=output, sep='\t')
            print(file=output)

        if counts:
            for kmer, count in zip(kmerset.kmers, kmerset.counts):
                print(kmertools.kmer_string(kmerset.k, int(kmer)), count,
                      sep='\t', file=output)

        if histogram:
            kmerset.write_histogram(output)
            print(file=output)

        if entropy:
            print("Entropy", round(kmerset.entropy(), 2),
                  file=output, sep='\t')


class PlotSubcommand(Subcommand):
    """
    Generate plots for a given k-mer set.
    """

    PLOT_TYPES = ('spectrum', )

    def register_arguments(self, subparser: argparse.ArgumentParser):
        subparser.add_argument(
            'kmerset',
            help="The k-mer set to load"
        )
        subparser.add_argument(
            '-o', '--output',
            help="Output filename (PNG preferred)."
        )
        subparser.add_argument(
            '-t', '--plot-type', choices=self.PLOT_TYPES,
            help="The kind of plot to generate."
        )

    def __call__(self, kmerset, output, plot_type, **kwargs):
        logger.info("Loading k-merset %s", kmerset)
        kmerset = kmertools.kmerset_from_hdf5(kmerset)

        if plot_type == 'spectrum':
            thresholds = kmerset.spectrum_min_max()
            if thresholds:
                kmerset.plot_spectrum(output, thresholds[2])
            else:
                kmerset.plot_spectrum(output)

            logger.info("Created k-mer spectrum plot in file %s", output)


class KmerizeSubcommand(Subcommand):
    """K-merize a given reference sequence or a sample read dataset."""

    def register_arguments(self, subparser: argparse.ArgumentParser):
        subparser.add_argument(
            'sequences', nargs='+',
            help='Input sequence files (fasta or fastq by default; optionally '
                 'compressed with gz or bz2)')
        subparser.add_argument(
            "-k", "--k", type=int, default=23,
            help="K-mer size (default %(default)s)",
        )
        subparser.add_argument(
            "-o", "--output",
            help="Filename of the output HDF5."
        )
        subparser.add_argument(
            "-f", "--fingerprint", action="store_true",
            help="Compute and save min-hash fingerprint (sketch)."
        )
        subparser.add_argument(
            '-s', '--sketch-fraction', type=float, default=0.01,
            help="Fraction of k-mers to keep for a minhash sketch. Default: "
                 "%(default)s"
        )
        subparser.add_argument(
            "-F", "--filter", action="store_true",
            help="Filter output kmers based on kmer spectrum (to prune "
                 "sequencing errors)"
        )
        subparser.add_argument(
            "-l", "--limit",
            help="Only process about this many kmers (can have suffix of M or"
                 " G)"
        )
        subparser.add_argument(
            "-p", "--prune",
            help="Prune singletons after accumulating this (can have suffix "
                 "of M or G)"
        )

    def __call__(self, k, sequences, output, limit=None, prune=None,
                 fingerprint=False, sketch_fraction=0.002, filter=False,
                 **kwargs):

        kmerset = kmertools.KmerSet(k)

        limit = utils.parse_num_suffix(limit)
        prune = utils.parse_num_suffix(prune)

        for seq in sequences:
            logger.info('K-merizing file %s...', seq)
            kmerset.kmerize_file(seq, limit=limit, prune=prune)

        if filter:
            thresholds = kmerset.spectrum_filter()
            if thresholds:
                logger.info("Filtered kmerset. Only k-mers within frequency "
                            "range [%d, %d] are kept.", *thresholds)

        if fingerprint:
            kmerset.min_hash(sketch_fraction)

        logger.info("Writing k-merset to %s", output)
        kmerset.save(output, compress=True)


class KmersimSubCommand(Subcommand):
    """
    Compare k-mer sets with each other. Both all-vs-all and one-vs-all is
    supported.
    """

    def register_arguments(self, subparser: argparse.ArgumentParser):
        compare_group = subparser.add_mutually_exclusive_group(required=True)
        compare_group.add_argument(
            '-a', '--all-vs-all', action="store_true", default=False,
            help="Perform all-vs-all comparisons for the given k-mer sets. "
                 "Either --all-vs-all is required or --sample."
        )
        compare_group.add_argument(
            '-s', '--sample', metavar='FILE',
            help="Perform one-vs-all comparisons with the given filename as "
                 "sample. Either --all-vs-all is required or --sample."
        )

        subparser.add_argument(
            '-f', '--fingerprint', action="store_true", default=False,
            required=False,
            help="Use min-hash fingerprint instead of full k-mer set."
        )
        subparser.add_argument(
            '-S', '--scoring', choices=list(comparison.SCORING_METHODS.keys()),
            default="jaccard", required=False,
            help="The scoring metric to use (default: jaccard)."
        )
        subparser.add_argument(
            '-F', '--fraction', action="store_true", default=False,
            required=False,
            help="Output numerator and denominator separately instead of "
                 "evaluating the division."
        )
        subparser.add_argument(
            '-t', '--threads', type=int, default=1, required=False,
            help="Use multiple processes the compute the similarity scores ("
                 "default 1)."
        )
        subparser.add_argument(
            '-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
            metavar='FILE',
            help="File to write the results (default: standard output)."
        )
        subparser.add_argument(
            'strains', nargs='+',
            help="Filenames of k-mer set HDF5 files."
        )

    def _do_compare(self, sets, scoring):
        try:
            set1, set2 = sets

            name1, data1 = set1
            name2, data2 = set2

            logging.info("Comparing %s vs %s...", name1, name2)

            similarity = comparison.similarity_score(data1, data2, scoring)

            return name1, name2, similarity
        except KeyboardInterrupt:
            pass

    def __call__(self, strains, output, all_vs_all=False, sample=None,
                 fingerprint=False, scoring="jaccard", threads=1,
                 fraction=False, **kwargs):

        # First, load all K-mer sets
        loader = (kmertools.load_fingerprint if fingerprint
                  else kmertools.load_kmers)
        logger.info("Loading %d k-mer sets...", len(strains))
        kmer_data = [(kmertools.name_from_path(strain), loader(strain))
                     for strain in strains]
        logging.info("Done.")

        to_compute_iter = None
        if sample:
            sample_data = (kmertools.name_from_path(sample), loader(sample))
            logger.info("Start %s vs all comparison...", sample_data[0])

            to_compute_iter = (
                (sample_data, strain_data) for strain_data in kmer_data
            )
        elif all_vs_all:
            if scoring == "reference":
                raise ValueError("'reference' scoring metric is meaningless in"
                                 " all-vs-all mode.")

            logger.info("Start computing pairwise similarities...")
            to_compute_iter = itertools.combinations(kmer_data, 2)

        if threads > 1:
            pool = multiprocessing.Pool(threads)

            scores = list(pool.imap_unordered(
                functools.partial(self._do_compare, scoring=scoring),
                to_compute_iter,
                chunksize=2**5
            ))
        else:
            scores = list(map(
                functools.partial(self._do_compare, scoring=scoring),
                to_compute_iter
            ))

        logger.info("Done.")

        # Sort results
        scores = sorted(scores, key=lambda e: e[2][0] / e[2][1], reverse=True)

        # Write results
        logger.info("Writing results...")
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        for name1, name2, (numerator, denominator) in scores:
            if fraction:
                writer.writerow((name1, name2, numerator, denominator,
                                 "{:.5f}".format(numerator / denominator)))
            else:
                writer.writerow((name1, name2, "{:.5f}".format(
                    numerator / denominator)))

        logger.info("Done.")


class ClusterSubcommand(Subcommand):
    """
    Group k-mer sets that are very similar to each other together.
    """

    def register_arguments(self, subparser: argparse.ArgumentParser):
        subparser.add_argument(
            '-c', '--cutoff', type=float, default=0.95,
            help="Minimum similarity between two sets to group them together."
        )
        subparser.add_argument(
            '-i', '--similarity-scores', type=argparse.FileType('r'),
            default=sys.stdin, metavar='FILE',
            help="The file with the similarity scores between sets (the "
                 "output of 'strainge compare --all-vs-all'). Defaults to"
                 " standard input."
        )
        subparser.add_argument(
            '-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
            metavar='FILE',
            help="The file where the list of sets to keep after clustering "
                 "gets written. Defaults to standard output."
        )
        subparser.add_argument(
            '--clusters-out', type=argparse.FileType('w'), default=None,
            required=False, metavar='FILE',
            help="Output an optional tab separated file with all clusters and "
                 "their entries."
        )
        subparser.add_argument(
            'kmersets', nargs='+', metavar='kmerset',
            help="The list of HDF5 filenames of k-mer sets to cluster."
        )

    def __call__(self, kmersets, similarity_scores, output, cutoff=0.95,
                 clusters_out=None, **kwargs):

        clusters = {}
        clustered_to = {}
        directories = {kmertools.name_from_path(f): os.path.dirname(f) for f in
                       kmersets}
        keep = set(directories.keys())

        num_clusters = 0
        logger.info("Start reading similarity scores...")
        for lineno, line in enumerate(similarity_scores):
            parts = line.strip().split()
            if len(parts) != 3:
                logger.warning("Line %d - invalid format: contains %d "
                               "elements instead of the expected 3 elements, "
                               "skipping.",
                               lineno, len(parts))
                continue

            similarity = float(parts[2])
            if similarity < cutoff:
                # We expect the list of similarity scores to be ordered,
                # when the similarity drops below the cutoff, we can stop
                logger.info("Similarity score dropped below cutoff (%.5f<%g)",
                            similarity, cutoff)
                break

            set1 = parts[0]
            set2 = parts[1]

            logger.info("%s is %.2f%% similar to %s", set1, similarity*100,
                        set2)

            if set1 in clustered_to and set2 not in clustered_to:
                cluster_id = clustered_to[set1]
                clusters[cluster_id].add(set2)
                clustered_to[set2] = cluster_id

                logger.debug("Added %s to existing cluster %d", set2,
                             cluster_id)
            elif set1 not in clustered_to and set2 in clustered_to:
                cluster_id = clustered_to[set2]
                clusters[cluster_id].add(set1)
                clustered_to[set1] = cluster_id

                logger.debug("Added %s to existing cluster %d", set1,
                             cluster_id)
            elif set1 in clustered_to and set2 in clustered_to:
                cid1 = clustered_to[set1]
                cid2 = clustered_to[set2]

                if cid1 != cid2:
                    # Merge the smaller cluster to the larger cluster
                    if len(clusters[cid1]) >= len(clusters[cid2]):
                        clusters[cid1].update(clusters[cid2])
                        for s in clusters[cid2]:
                            clustered_to[s] = cid1

                        del clusters[cid2]
                        logger.debug("Merged cluster %d into %d", cid2, cid1)
                    else:
                        clusters[cid2].update(clusters[cid1])
                        for s in clusters[cid1]:
                            clustered_to[s] = cid2

                        del clusters[cid1]
                        logger.debug("Merged cluster %d into %d", cid1, cid2)
            else:
                # Create new cluster
                clusters[num_clusters] = {set1, set2}
                clustered_to[set1] = num_clusters
                clustered_to[set2] = num_clusters
                logger.debug("Created new cluster with ID %d", num_clusters)

                num_clusters += 1

        logger.info("Clustering done.")

        writer = None
        if clusters_out:
            writer = csv.writer(clusters_out, delimiter="\t",
                                lineterminator="\n")

        logger.info("Figure out which sets to keep (by sorting the entries by "
                    "number of scaffolds in the FASTA)...")

        def _kmerset_num_scaffolds(kmerset):
            filename = os.path.join(directories[kmerset], kmerset) + ".hdf5"
            with h5py.File(filename, 'r') as h5:
                if 'nSeqs' not in h5.attrs:
                    raise KeyError(
                        f"k-mer set '{filename}' does not contain number of "
                        "sequences in its HDF5 file. Key 'nSeqs' missing in "
                        "HDF5 attributes."
                    )

                return int(h5.attrs['nSeqs'])

        for cluster in sorted(clusters):
            sets = clusters[cluster]
            sorted_sets = sorted(sets, key=_kmerset_num_scaffolds)

            if clusters_out:
                items = [cluster]
                items.extend(sorted_sets)
                writer.writerow(items)

            sorted_sets.pop(0)
            keep.difference_update(sorted_sets)

        logger.info("Done. After clustering %d/%d genomes remain. Writing "
                    "output...", len(keep), len(kmersets))
        for kmerset in keep:
            print(os.path.join(directories[kmerset], kmerset) + ".hdf5",
                  file=output)


class CreateDBSubcommand(Subcommand):
    """
    Create pan-genome database in HDF5 format from a list of k-merized
    strains.
    """

    def register_arguments(self, subparser: argparse.ArgumentParser):
        subparser.add_argument(
            '-o', '--output',
            help="Pan-genome database output HDF5 file."
        )
        subparser.add_argument(
            '-F', '--fingerprint', action="store_true", default=False,
            help="Create fingerprint for the pan-genome database."
        )
        subparser.add_argument(
            '-f', '--from-file', type=argparse.FileType('r'), default=None,
            metavar='FILE',
            help="Read list of HDF5 filenames to include in the database from "
                 "a given file (use '-' to denote standard input). This is in "
                 "addition to any k-merset given as positional argument."
        )
        subparser.add_argument(
            'kmersets', metavar='kmerset', nargs='*',
            help="The HDF5 filenames of the kmerized reference strains."
        )

    def __call__(self, kmersets, from_file, output, fingerprint=False,
                 **kwargs):
        if from_file:
            for line in from_file:
                kmersets.append(line.strip())

        if not kmersets:
            logger.error("No k-mer sets given and nothing read from the given "
                         "file, stopping.")
            return 1

        pankmerset = None
        with h5py.File(output, 'w') as h5:
            for fname in kmersets:
                name = kmertools.name_from_path(fname)
                kset = kmertools.kmerset_from_hdf5(fname)
                logger.info("Adding k-merset %s", name)

                strain_group = h5.create_group(name)
                kset.save_hdf5(strain_group, compress="gzip")

                if not pankmerset:
                    pankmerset = kset
                else:
                    pankmerset = pankmerset.merge_kmerset(kset)

            if fingerprint:
                pankmerset.min_hash()

            logger.info("Saving pan-genome database")
            pankmerset.save_hdf5(h5, compress="gzip")
            logger.info("Done.")
