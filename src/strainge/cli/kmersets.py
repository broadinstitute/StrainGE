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

import sys
import csv
import logging
import argparse
import functools
import itertools
import multiprocessing
from pathlib import Path

import h5py
import numpy
import pandas

from strainge.cli.registry import Subcommand
from strainge import kmertools, utils, comparison

from strainge import cluster

logger = logging.getLogger()


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
            "-o", "--output", required=True,
            help="Filename of the output HDF5."
        )
        subparser.add_argument(
            '-f', '--fingerprint-fraction', type=float, default=kmertools.DEFAULT_FINGERPRINT_FRACTION,
            help="Fraction of k-mers to keep for a minhash sketch. Default: "
                 "%(default)s. No fingerprint will be created if set to zero."
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
                 fingerprint_fraction=kmertools.DEFAULT_FINGERPRINT_FRACTION, filter=False,
                 **kwargs):

        kmerset = kmertools.KmerSet(k)

        limit = utils.parse_num_suffix(limit)
        prune = utils.parse_num_suffix(prune)

        if not output:
            logger.error("No output filename given! Please specify the output file with `-o`.")
            return 1

        for seq in sequences:
            logger.info('K-merizing file %s...', seq)
            kmerset.kmerize_file(seq, limit=limit, prune=prune)

        if filter:
            thresholds = kmerset.spectrum_filter()
            if thresholds:
                logger.info("Filtered kmerset. Only k-mers within frequency "
                            "range [%d, %d] are kept (mode %d).", thresholds[0], thresholds[2], thresholds[1])
            else:
                logger.warning("Unable to compute filter thresholds! Continuing without filtering...")

        if fingerprint_fraction:
            kmerset.min_hash(fingerprint_fraction)

        logger.info("Writing k-merset to %s", output)
        kmerset.save(output, compress=True)


class KmermergeSubcommand(Subcommand):
    """Merge k-mer set files."""

    def register_arguments(self, subparser: argparse.ArgumentParser):
        subparser.add_argument(
            'kmerfiles', nargs='+',
            help='Input KmerSet files to be merged (hdf5 files)'
        )
        subparser.add_argument(
            "-k", "--k", type=int, default=23,
            help="K-mer size (default %(default)s)",
        )
        subparser.add_argument(
            "-o", "--output",
            help="Filename of the output HDF5."
        )
        subparser.add_argument(
            '-f', '--fingerprint-fraction', type=float, default=kmertools.DEFAULT_FINGERPRINT_FRACTION,
            help="Fraction of k-mers to keep for a minhash sketch. Default: "
                 "%(default)s. No fingerprint will be created if set to zero."
        )

    def __call__(self, k, kmerfiles, output, limit=None, prune=None,
                 fingerprint_fraction=kmertools.DEFAULT_FINGERPRINT_FRACTION, filter=False,
                 **kwargs):

        kmerset = None

        for ksfile in kmerfiles:
            logger.info('Merging KmerSet file %s...', ksfile)
            ks = kmertools.kmerset_from_hdf5(ksfile)
            assert ks.k == k, "Incompatible kmer size {}".format(ks.k)
            kmerset = kmerset.merge_kmerset(ks) if kmerset else ks

        if fingerprint_fraction:
            kmerset.min_hash(fingerprint_fraction)

        logger.info("Writing k-merset to %s", output)
        kmerset.save(output, compress=True)


class KmersimRunner:
    """
    This class helps running the comparisons using multiprocessing. We cache
    kmersets to make sure we don't have to reload them each time. It's
    possible, however, that different processes have to cache the same
    kmerset. By using a large enough chunksize (see
    `KmersimSubcommand.__call__`) we aim to prevent that from happening too
    often.
    """

    def __init__(self, kmersets, scoring, fingerprint):
        logger.info("Loading k-mer sets...")
        with h5py.File(kmersets[0], 'r') as h5:
            self.k = h5.attrs['k']

        load_func = kmertools.load_fingerprint if fingerprint else kmertools.load_kmers
        self.kmersets = {
            kmerset: load_func(kmerset, expect_k=self.k) for kmerset in kmersets
        }

        self.scoring = scoring

    def __call__(self, kmersets):
        try:
            set1, set2 = kmersets

            name1 = kmertools.name_from_path(set1)
            name2 = kmertools.name_from_path(set2)

            data1 = self.kmersets[set1]
            data2 = self.kmersets[set2]

            logger.info("Comparing %s vs %s...", name1, name2)

            scores = {
                metric: comparison.similarity_score(data1, data2, metric)
                for metric in self.scoring if metric != 'subset'
            }

            if 'jaccard' in scores:
                scores['ani'] = comparison.ani(self.k, scores['jaccard'])

            # Subset scores are not symmetric
            if 'subset' in self.scoring:
                scores['subset1'] = comparison.subset(data1, data2)
                scores['subset2'] = comparison.subset(data2, data1)

            return [name1, name2, scores]
        except KeyboardInterrupt:
            pass


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
            '-f', '--full-db', action="store_true",
            help="Use full k-mer set instead of min-hash fingerprint."
        )
        subparser.add_argument(
            '-S', '--scoring', choices=list(comparison.SCORING_METHODS.keys()),
            default=[], required=False, action="append",
            help="The scoring metric to use (default: jaccard). Can be used "
                 "multiple times to include multiple scoring metrics. Choices:"
                 " %(choices)s."
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

    def __call__(self, strains, output, all_vs_all=False, sample=None,
                 full_db=False, scoring=None, threads=1,
                 fraction=False, **kwargs):

        if not scoring:
            scoring = ["jaccard", "subset"]

        fingerprint = not full_db

        if sample:
            sample_name = kmertools.name_from_path(sample)
            logger.info("Start %s vs all comparison...", sample_name)

            to_compute_iter = (
                (sample, strain) for strain in strains
            )
        elif all_vs_all:
            if scoring == "reference":
                raise ValueError("'reference' scoring metric is meaningless in"
                                 " all-vs-all mode.")

            logger.info("Start computing pairwise similarities...")
            to_compute_iter = itertools.combinations(strains, 2)
        else:
            logger.error("Either --sample or --all-vs-all required.")
            return 1

        runner = KmersimRunner(strains, scoring, fingerprint)
        if threads > 1:
            pool = multiprocessing.Pool(threads)

            scores = list(pool.imap_unordered(runner, to_compute_iter,
                                              chunksize=2**16))
        else:
            scores = list(map(runner, to_compute_iter))

        logger.info("Done.")

        # Sort results
        # We use the first given scoring metric for sorting
        first_metric = scoring[0]
        scores = sorted(scores, key=lambda e: e[2][first_metric], reverse=True)

        # Write results
        logger.info("Writing results...")
        scoring_keys = scores[0][2].keys()
        fieldnames = ["kmerset1", "kmerset2", *scoring_keys]
        writer = csv.DictWriter(output, fieldnames, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for name1, name2, pair_score in scores:
            data = {metric: f"{value:.5f}"
                    for metric, value in pair_score.items()}
            data['kmerset1'] = name1
            data['kmerset2'] = name2

            writer.writerow(data)

        logger.info("Done.")


class ClusterSubcommand(Subcommand):
    """
    Group k-mer sets that are very similar to each other together.
    """

    def register_arguments(self, subparser: argparse.ArgumentParser):
        subparser.add_argument(
            '-c', '--cutoff', type=float, default=0.90,
            help="Minimum similarity between two sets to group them together."
        )
        subparser.add_argument(
            '-i', '--similarity-scores', type=argparse.FileType('r'),
            default=sys.stdin, metavar='FILE',
            help="The file with the similarity scores between kmersets (the "
                 "output of 'strainge compare --all-vs-all'). Defaults to"
                 " standard input."
        )
        subparser.add_argument(
            '-d', '--discard-contained', action="store_true",
            help="Discard k-mersets that are a subset of another k-merset. "
                 "Requires 'subset' scoring metric in the similarity scores "
                 "TSV files."
        )
        subparser.add_argument(
            '-C', '--contained-cutoff', type=float, default=0.99,
            help="Minimum fraction of kmers to be present in another genome "
                 "to discard it."
        )
        subparser.add_argument(
            '-w', '--warn-too-distant', type=float, default=85, metavar="ANI",
            help="Warn when including references that that seem too distantly related, which could indicate a "
                 "mislabeled reference genome. Default: %(default)s%% ANI."
        )
        subparser.add_argument(
            '-p', '--priorities', type=argparse.FileType('r'), default=None,
            metavar="FILE",
            help="An optional TSV file where the first column represents the "
                 "ID of a reference kmerset, and the second an integer "
                 "indicating the priority for clustering. References with "
                 "higher priority get precedence over references with lower "
                 "priority in the same cluster."
        )
        subparser.add_argument(
            '-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
            metavar='FILE',
            help="The file where the list of kmersets to keep after "
                 "clustering gets written. Defaults to standard output."
        )
        subparser.add_argument(
            '--clusters-out', type=argparse.FileType('w'), default=None,
            required=False, metavar='FILE',
            help="Output an optional tab separated file with all clusters and "
                 "their entries."
        )
        subparser.add_argument(
            'kmersets', nargs='+', metavar='kmerset', type=Path,
            help="The list of HDF5 filenames of k-mer sets to cluster."
        )

    def __call__(self, kmersets, similarity_scores, output,
                 discard_contained=False, priorities=None, cutoff=0.95, warn_too_distant=85,
                 contained_cutoff=0.99, clusters_out=None, **kwargs):
        label_to_path = {
            kset.stem: kset for kset in kmersets
        }
        labels = list(label_to_path.keys())
        label_ix = {label: i for i, label in enumerate(labels)}

        logger.info("Reading pairwise similarities...")
        similarities = pandas.read_csv(similarity_scores, sep='\t', comment='#').set_index(['kmerset1', 'kmerset2'])

        # Check for references too distant
        if 'ani' in similarities:
            sim_matrix = cluster.similarities_to_matrix(similarities, labels, 'ani')

            threshold = warn_too_distant / 100
            mean_sim = numpy.nanmean(sim_matrix.values, axis=1)

            ix = mean_sim < threshold
            too_distant = sim_matrix[mean_sim < threshold].index

            for ref, dist in zip(too_distant, mean_sim[ix]):
                logger.warning(f"%s average ANI to other references is {dist*100:.1f}%%, possibly mislabeled "
                               "reference?", ref)

        logger.info("Clustering genomes...")
        clusters = cluster.cluster_genomes(similarities, labels, cutoff)

        ref_priorities = {}
        if priorities:
            reader = csv.reader(priorities, delimiter='\t')
            ref_priorities = {
                row[0].strip(): int(row[1].strip()) for row in reader if
                len(row) == 2
            }

        logger.info("Picking a representive genome per cluster...")
        count = 0
        to_keep = []
        for ix, sorted_entries in cluster.pick_representative(
                clusters, similarities, ref_priorities):
            to_keep.append(sorted_entries[0])

            if clusters_out:
                print(*sorted_entries, sep='\t', file=clusters_out)

            count += 1

        similarities = similarities.loc[(to_keep, to_keep), :].copy()

        exclude = set()
        if discard_contained and 'subset1' not in similarities.columns:
            logger.warning("No 'subset' score in similarities file. Can't "
                           "discard k-mersets that are subsets of another. "
                           "Run `strainge kmersim` with both '--scoring "
                           "jaccard' and '--scoring subset'")
        elif discard_contained:
            subset_matrix = cluster.similarities_to_matrix(similarities, to_keep, 'subset')
            subset_max = numpy.nanmax(subset_matrix, axis=1)

            exclude = subset_matrix[subset_max >= contained_cutoff].index
            logger.info("Excluding %d genomes because they're contained for at least "
                        "%d%% in another genome.", len(exclude),
                        contained_cutoff * 100)

            to_keep = [l for l in to_keep if l not in exclude]

        for l in to_keep:
            print(label_to_path[l], file=output)

        logger.info("Done. After clustering %d/%d genomes remain.",
                    count, len(kmersets))


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

    def __call__(self, kmersets, from_file, output,
                 **kwargs):
        if from_file:
            for line in from_file:
                kmersets.append(line.strip())

        if not kmersets:
            logger.error("No k-mer sets given and nothing read from the given "
                         "file, stopping.")
            return 1

        pankmerset = None
        fingerprint = True
        panfingerprint = None
        fpf = None

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
                    if kset.fingerprint is not None:
                        fp = kset.fingerprint_as_kmerset()
                        if not panfingerprint:
                            panfingerprint = fp
                            fpf = fp.fingerprint_fraction
                        else:
                            panfingerprint = panfingerprint.merge_kmerset(fp)
                    else:
                        logger.warning("Not all input kmersets have fingerprints, so no pan genome fingerprint will be generated")
                        fingerprint = False

            if fingerprint and panfingerprint is not None:
                logger.info("Adding pan-genome fingerprint (%d distinct kmers)", panfingerprint.kmers.size)
                pankmerset.fingerprint = panfingerprint.kmers
                pankmerset.fingerprint_counts = panfingerprint.counts
                pankmerset.fingerprint_fraction = fpf

            logger.info("Saving pan-genome database")
            pankmerset.save_hdf5(h5, compress="gzip")
            logger.info("Done.")
