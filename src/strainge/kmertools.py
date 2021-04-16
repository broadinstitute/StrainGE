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
import logging

import h5py
import pysam
import skbio
import numpy as np
import matplotlib.pyplot as plt

from strainge import kmerizer
from strainge.io.utils import open_compressed, read_fastq

logger = logging.getLogger(__name__)

DEFAULT_K = 23
DEFAULT_FINGERPRINT_FRACTION = 0.01
OLD_FINGERPRINT_FRACTION = 0.002

A = 0
C = 1
G = 2
T = 3

BASES = "ACGT"


def kmer_string(k, kmer):
    seq = ''.join([BASES[(kmer >> (2 * k)) & 3] for k in range(k - 1, -1, -1)])
    return seq


def iter_sequences_bam(bamfile):
    """Iterate over sequences in a BAM file. Only outputs the sequence, useful
    for kmerizing."""

    bam = pysam.AlignmentFile(bamfile, check_header=False, check_sq=False)

    seq_iter = iter(bam.fetch(until_eof=True))
    seq_iter = filter(lambda r: not r.is_qcfail, seq_iter)

    yield from (seq.seq.encode('utf-8') for seq in seq_iter)

    bam.close()


def iter_sequences_fasta(f):
    """Use scikit-bio to iterate over FASTA sequences."""

    yield from (str(seq) for seq in skbio.io.read(f, "fasta"))


def iter_sequences_fastq(f):
    """Use Heng Li's fast FASTQ reader to iterate over reads"""

    yield from (r[1] for r in read_fastq(f))


def open_seq_file(file_name):
    """
    Iterate over sequences present in either a BAM file, FASTA file, or FASTQ
    file.

    Assumes fasta unless ".fastq" or ".fq" in the file name.

    Parameters
    ----------
    file_name : str
        The file to open

    Yields
    ------
    str
        Each sequence present in the given file
    """

    components = file_name.split('.')

    if "bam" in components:
        yield from iter_sequences_bam(file_name)
    else:
        with open_compressed(file_name) as f:
            if "fastq" in components or "fq" in components:
                yield from iter_sequences_fastq(f)
            else:
                yield from iter_sequences_fasta(f)


def load_hdf5(file_path, thing, expect_k=None):
    with h5py.File(file_path, 'r') as h5:
        hdf5_type = h5.attrs['type']

        if isinstance(hdf5_type, bytes):
            hdf5_type = hdf5_type.decode()

        if hdf5_type != "KmerSet":
            raise ValueError("The HDF5 file is not a KmerSet, unexpected type:"
                             " '{}'".format(h5.attrs['type']))

        k = h5.attrs['k']
        if expect_k is not None and expect_k != k:
            raise ValueError(f"The loaded kmerset has not the expected k-mer size! Expected: {expect_k}, actual: {k}")

        return np.array(h5[thing])


def load_kmers(file_name, expect_k=None):
    return load_hdf5(file_name, "kmers", expect_k)


def load_counts(file_name, expect_k=None):
    return load_hdf5(file_name, "counts", expect_k)


def load_fingerprint(file_name, expect_k=None):
    return load_hdf5(file_name, "fingerprint", expect_k)


def name_from_path(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]


def kmerset_from_hdf5(file_path):
    if not file_path.endswith(".hdf5"):
        file_path += ".hdf5"

    with h5py.File(file_path, 'r') as h5:
        hdf5_type = h5.attrs['type']

        if isinstance(hdf5_type, bytes):
            hdf5_type = hdf5_type.decode()

        assert hdf5_type == "KmerSet", "Not a KmerSet file!"
        kset = KmerSet(h5.attrs['k'])

        if "fingerprint_fraction" in h5.attrs:
            kset.fingerprint_fraction = h5.attrs["fingerprint_fraction"]
        if "fingerprint" in h5:
            kset.fingerprint = np.array(h5["fingerprint"])
            if not kset.fingerprint_fraction:
                kset.fingerprint_fraction = OLD_FINGERPRINT_FRACTION
        if "fingerprint_counts" in h5:
            kset.fingerprint_counts = np.array(h5["fingerprint_counts"])

        if "kmers" in h5:
            kset.kmers = np.array(h5["kmers"])
        if "counts" in h5:
            kset.counts = np.array(h5["counts"])

    return kset


def kmerset_from_file(file_path, k=DEFAULT_K):
    return kmerset_from_hdf5(file_path)


def similarity_score(kmers1, kmers2, scoring="jaccard"):
    """Compute Jaccard similarity index"""
    # count of kmers in common
    intersection = float(kmerizer.count_common(kmers1, kmers2))
    if scoring == "jaccard":
        # Use Jaccard similarity index
        score = intersection / (kmers1.size + kmers2.size - intersection)
    elif scoring == "minsize":
        # Use intersection / min_size (proper subset scores 1.0)
        score = intersection / min(kmers1.size, kmers2.size)
    elif scoring == "meansize":
        # Use mean size in denominator (used in Mash)
        score = intersection / ((kmers1.size + kmers2.size) / 2)
    elif scoring == "maxsize":
        # Use intersection / max_size (proper subset scores min/max)
        score = intersection / max(kmers1.size, kmers2.size)
    elif scoring == "reference":
        # Use intersection / size of reference (useful for comparing reads to
        # assembled references)
        score = intersection / kmers2.size
    else:
        assert scoring in (
            "jaccard", "minsize", "maxsize", "meansize", "reference"), \
            "unknown scoring method"
    return score


def similarity_numerator_denominator(kmers1, kmers2, scoring="jaccard"):
    """Compute Jaccard similarity index"""
    # count of kmers in common
    intersection = float(kmerizer.count_common(kmers1, kmers2))
    if scoring == "jaccard":
        # Use Jaccard similarity index
        denom = (kmers1.size + kmers2.size - intersection)
    elif scoring == "minsize":
        # Use intersection / min_size (proper subset scores 1.0)
        denom = min(kmers1.size, kmers2.size)
    elif scoring == "maxsize":
        # Use intersection / max_size (proper subset scores min/max)
        denom = max(kmers1.size, kmers2.size)
    elif scoring == "reference":
        # Use intersection / size of reference (useful for comparing reads to
        # assembled references)
        denom = kmers2.size
    else:
        assert scoring in ("jaccard", "minsize", "maxsize"), \
            "unknown scoring method"
    return intersection, denom


def build_kmer_count_matrix(kmersets):
    """Build a big matrix with kmer counts from a list of kmersets.

    Each column will represent a single k-mer set and each row a k-mer. This
    will effectively merge all kmersets to a single matrix.

    Parameters
    ----------
    kmersets : List[KmerSet]
        List of `KmerSet` objects to build the matrix from.

    Returns
    -------
    Tuple[List[kmer_t], array]
        This function returns a tuple with two elements: the first element is
        a list of k-mers, i.e. the labels for the rows of the matrix, and the
        second element is the matrix itself.
    """

    # Defer to our C++ extension
    return kmerizer.build_kmer_count_matrix([
        (kmerset.kmers, kmerset.counts) for kmerset in kmersets
    ])


class KmerSet(object):
    """
    Holds array of kmers and their associated counts & stats.
    """

    def __init__(self, k=DEFAULT_K):
        self.k = k
        # data arrays
        self.kmers = None
        self.counts = None
        self.fingerprint = None
        self.fingerprint_counts = None
        self.fingerprint_fraction = None

        self.singletons = None

        # stats from kmerizing, if appropriate
        self.n_seqs = 0
        self.n_bases = 0
        self.n_kmers = 0

    def __eq__(self, other):
        return (self.k == other.k
                and np.array_equal(self.fingerprint, other.fingerprint)
                and np.array_equal(self.kmers, other.kmers)
                and np.array_equal(self.counts, other.counts))

    def kmerize_file(self, file_name, batch_size=100000000, verbose=True,
                     limit=0, prune=0):
        seq_file = open_seq_file(file_name)
        batch = np.empty(batch_size, dtype=np.uint64)

        n_seqs = 0
        n_bases = 0
        n_kmers = 0
        pruned = False

        for seq in seq_file:
            n_seqs += 1
            seq_length = len(seq)
            n_bases += seq_length
            if n_kmers + seq_length > batch_size:
                self.process_batch(batch, n_seqs, n_bases, n_kmers, verbose)
                if limit and self.n_kmers > limit:
                    break
                if prune and self.singletons > prune:
                    self.prune_singletons(verbose)
                    pruned = True
                n_seqs = 0
                n_bases = 0
                n_kmers = 0

            n_kmers += kmerizer.kmerize_into_array(self.k, seq, batch, n_kmers)
            if limit and self.n_kmers + n_kmers >= limit:
                break

        self.process_batch(batch, n_seqs, n_bases, n_kmers, verbose)
        if pruned:
            self.prune_singletons(verbose)

    def kmerize_seq(self, seq):
        kmers = kmerizer.kmerize(self.k, seq)
        self.n_seqs += 1
        self.n_bases += len(seq)
        self.n_kmers = kmers.size
        self.kmers, self.counts = np.unique(kmers, return_counts=True)

    def process_batch(self, batch, nseqs, nbases, nkmers, verbose):
        self.n_seqs += nseqs
        self.n_bases += nbases
        self.n_kmers += nkmers

        new_kmers, new_counts = np.unique(batch[:nkmers], return_counts=True)

        if self.kmers is None:
            self.kmers = new_kmers
            self.counts = new_counts
        else:
            self.kmers, self.counts = kmerizer.merge_counts(
                self.kmers, self.counts, new_kmers, new_counts)

        self.singletons = np.count_nonzero(self.counts == 1)
        if verbose:
            self.print_stats()

    def prune_singletons(self, verbose=False):
        keepers = self.counts > 1
        self.kmers = self.kmers[keepers]
        self.counts = self.counts[keepers]
        logger.debug("Pruned singletons: %d distinct k-mers remain",
                     self.kmers.size)

    def merge_kmerset(self, other):
        """Create new KmerSet by merging this with another"""
        new_set = KmerSet(self.k)
        new_set.kmers, new_set.counts = kmerizer.merge_counts(
            self.kmers, self.counts, other.kmers, other.counts)
        return new_set

    def intersect(self, kmers):
        """
        Compute intersection with given kmers
        :param kmers: kmers to keep
        :return: reduced version of self
        """

        ix = kmerizer.intersect_ix(self.kmers, kmers)
        self.counts = self.counts[ix]
        self.kmers = self.kmers[ix]

        return self

    def exclude(self, kmers):
        """
        Return this KmerSet with excluded kmers removed.
        :param kmers: kmers to exclude
        :return: reduced version of self
        """
        new_kmers = kmerizer.diff(self.kmers, kmers)

        ix = kmerizer.intersect_ix(self.kmers, new_kmers)
        self.counts = self.counts[ix]
        self.kmers = new_kmers

        return self

    def mutual_intersect(self, other):
        """
        Compute intersection of two kmer sets and reduce both to their common
        kmers. BOTH sets are modified!

        :param other: other KmerSet
        :return: reduced self
        """
        ix = kmerizer.intersect_ix(self.kmers, other.kmers)
        self.kmers = self.kmers[ix]
        self.counts = self.counts[ix]

        ix = kmerizer.intersect_ix(other.kmers, self.kmers)
        other.kmers = other.kmers[ix]
        other.counts = other.counts[ix]

        return self

    def print_stats(self):
        logger.info("Seqs %d, bases %d, kmers: %d, distinct: %d, singletons: "
                    "%d", self.n_seqs, self.n_bases, self.n_kmers,
                    self.kmers.size, self.singletons)

    def min_hash(self, frac=DEFAULT_FINGERPRINT_FRACTION):
        nkmers = int(round(self.kmers.size * frac))
        order = kmerizer.fnvhash_kmers(self.k, self.kmers).argsort()[:nkmers]
        self.fingerprint = self.kmers[order]
        self.fingerprint.sort()

        ix = kmerizer.intersect_ix(self.kmers, self.fingerprint)
        self.fingerprint_counts = self.counts[ix]

        self.fingerprint_fraction = frac

        return self.fingerprint


    def fingerprint_override(self):
        """Use the fingerprint instead of the full KmerSet"""
        assert self.fingerprint is not None, "Can't use fingerprint if it's not there!"
        self.kmers = self.fingerprint
        self.counts = self.fingerprint_counts


    def fingerprint_as_kmerset(self):
        assert self.fingerprint is not None

        kset = KmerSet(k=self.k)
        kset.kmers = self.fingerprint
        kset.fingerprint_fraction = self.fingerprint_fraction

        if self.fingerprint_counts is not None:
            kset.counts = self.fingerprint_counts
        else:
            kset.counts = np.ones_like(kset.kmers, dtype=np.uint64)

        return kset

    def freq_filter(self, min_freq=1, max_freq=None):
        condition = (self.counts >= min_freq)
        if max_freq:
            condition &= (self.counts <= max_freq)
        self.kmers = self.kmers[condition]
        self.counts = self.counts[condition]

    def spectrum(self):
        return np.unique(self.counts, return_counts=True)

    def spectrum_min_max(self, delta=.5, max_copy_number=20):
        freq, counts = self.spectrum()
        min_index = 0
        max_index = 0
        have_min = False
        have_max = False
        last_freq = 0
        for i in range(freq.size):
            count = counts[i]
            zero = freq[i] > last_freq + 1
            if have_max and (
                    zero or freq[i] > freq[max_index] * max_copy_number):
                break
            if have_min:
                if count > counts[max_index]:
                    max_index = i
                if count < counts[max_index] * (1 - delta):
                    have_max = True
            elif count > 1000 and count > counts[min_index] * (1 + delta):
                have_min = True
            elif zero or count < counts[min_index]:
                min_index = i
                max_index = i
            elif counts[i] < counts[min_index]:
                min_index = max_index = i
            last_freq = freq[i]

        if (min_index and max_index
                and counts[max_index] > counts[min_index] * (1 + delta)):
            return freq[min_index], freq[max_index], freq[i - 1]

        return None

    def spectrum_filter(self, max_copy_number=20):
        thresholds = self.spectrum_min_max()
        if thresholds:
            self.freq_filter(thresholds[0], thresholds[2])
        return thresholds

    def plot_spectrum(self, file_name=None, max_freq=None):
        # to get kmer profile, count the counts!
        spectrum = self.spectrum()
        plt.semilogy(spectrum[0], spectrum[1])
        plt.grid = True
        if max_freq:
            plt.xlim(0, max_freq)
        plt.xlabel("Kmer Frequency")
        plt.ylabel("Number of Kmers")
        if file_name:
            plt.savefig(file_name)
        else:
            plt.show()

    def write_histogram(self, file_obj):
        spectrum = self.spectrum()
        for i in range(spectrum[0].size):
            print("%d\t%d" % (spectrum[0][i], spectrum[1][i]), file=file_obj)

    def entropy(self):
        """Calculate Shannon entropy in bases"""
        if self.counts is None:
            return 0.0
        total = float(self.counts.sum())
        probs = self.counts / total
        return (-(probs * np.log2(probs)).sum()) / 2

    def save_hdf5(self, h5, compress=None):
        h5.attrs["type"] = "KmerSet"
        h5.attrs["k"] = self.k
        h5.attrs["nSeqs"] = self.n_seqs

        if self.fingerprint is not None:
            h5.create_dataset("fingerprint", data=self.fingerprint,
                              compression=compress)
        if self.fingerprint_counts is not None:
            h5.create_dataset("fingerprint_counts",
                              data=self.fingerprint_counts,
                              compression=compress)
        if self.fingerprint_fraction is not None:
            h5.attrs["fingerprint_fraction"] = self.fingerprint_fraction

        if self.kmers is not None:
            h5.create_dataset("kmers", data=self.kmers, compression=compress)
        if self.counts is not None:
            h5.create_dataset("counts", data=self.counts, compression=compress)

    def save(self, file_name, compress=None):
        """Save in HDF5 file format"""
        if compress is True:
            compress = "gzip"
        if not file_name.endswith(".hdf5"):
            file_name += ".hdf5"
        with h5py.File(file_name, 'w') as h5:
            self.save_hdf5(h5, compress)

    def load_hdf5(self, h5):
        h5_type = h5.attrs['type']

        # Support for HDF5 files generated in previous versions under Python 2
        if isinstance(h5_type, bytes):
            h5_type = h5_type.decode()

        if h5_type != "KmerSet":
            raise ValueError("The HDF5 file is not a KmerSet, unexpected type:"
                             " '{}'".format(h5_type))

        self.k = int(h5.attrs['k'])
        if 'nSeqs' in h5.attrs:
            self.n_seqs = int(h5.attrs['nSeqs'])

        if "fingerprint_fraction" in h5.attrs:
            self.fingerprint_fraction = h5.attrs["fingerprint_fraction"]
        if "fingerprint" in h5:
            self.fingerprint = np.array(h5["fingerprint"])
            if not self.fingerprint_fraction:
                self.fingerprint_fraction = OLD_FINGERPRINT_FRACTION
        if "fingerprint_counts" in h5:
            self.fingerprint_counts = np.array(h5["fingerprint_counts"])

        if "kmers" in h5:
            self.kmers = np.array(h5["kmers"])
        if "counts" in h5:
            self.counts = np.array(h5["counts"])

    def load(self, file_name):
        with h5py.File(file_name, 'r') as h5:
            self.load_hdf5(h5)

    def copy(self):
        new_kmerset = KmerSet(self.k)
        new_kmerset.kmers = self.kmers.copy()
        new_kmerset.counts = self.counts.copy()

        if self.fingerprint is not None:
            new_kmerset.fingerprint = self.fingerprint.copy()
            new_kmerset.fingerprint_counts = self.fingerprint_counts.copy()

        new_kmerset.fingerprint_fraction = self.fingerprint_fraction

        new_kmerset.singletons = self.singletons
        new_kmerset.n_seqs = self.n_seqs
        new_kmerset.n_kmers = self.n_kmers
        new_kmerset.n_bases = self.n_bases

        return new_kmerset
