import os
import gzip
import bz2
import h5py
import pysam
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import kmerizer

DEFAULT_K = 23

A = 0
C = 1
G = 2
T = 3

BASES = "ACGT"


def kmerString(k, kmer):
    seq = ''.join([BASES[(kmer >> (2 * k)) & 3] for k in xrange(k - 1, -1, -1)])
    #seqrc = Seq.reverse_complement(seq)
    #return ("%x" % kmer, seq, seqrc)
    return seq


def openSeqFile(fileName):
    """
    Open a sequence file with SeqIO; can be fasta or fastq with optional gz or bz2 compression.
    Assumes fasta unless ".fastq" or ".fq" in the file name.
    :param fileName:
    :return: SeqIO.parse object
    """

    components = fileName.split('.')

    if "bam" in components:
        file = pysam.AlignmentFile(fileName, "rb", check_header=False, check_sq=False)

        # generator for sequences in bam
        def bamSequences():
            for read in file.fetch(until_eof=True):
                if not read.is_qcfail:
                    yield read

        return bamSequences()

    if "bz2" in components:
        file = bz2.BZ2File(fileName, 'r')
    elif "gz" in components:
        file = gzip.GzipFile(fileName, 'r')
        SeqIO.parse(file, "fastq")
    else:
        file = open(fileName, 'r')
    if "fastq" in components or "fq" in components:
        fileType = "fastq"
    else:
        fileType = "fasta"
    return SeqIO.parse(file, fileType)

def loadNpz(fileName, thing):
    """
    :param fileName: A numpy npz file
    :param thing: A file element of the npz file to return
    :return: the the thing array
    """
    data = np.load(fileName)
    return data[thing]

def loadHdf5(filePath, thing):
    with h5py.File(filePath, 'r') as h5:
        assert h5.attrs["type"] == "KmerSet", "Not a KmerSet file!"
        return np.array(h5[thing])

def loadKmers(fileName):
    if fileName.endswith(".npz"):
        return loadNpz(fileName, "kmers")
    else:
        return loadHdf5(fileName, "kmers")

def loadCounts(fileName):
    if fileName.endswith(".npz"):
        return loadNpz(fileName, "counts")
    else:
        return loadHdf5(fileName, "counts")

def loadFingerprint(fileName):
    if fileName.endswith(".npz"):
        return loadNpz(fileName, "fingerprint")
    else:
        return loadHdf5(fileName, "fingerprint")

def nameFromPath(filePath):
    return os.path.splitext(os.path.basename(filePath))[0]


def kmerSetFromNpz(filePath, k):
    if not filePath.endswith(".npz"):
        filePath += ".npz"
    kset = KmerSet(k)
    data = np.load(filePath)
    if "fingerprint" in data:
        kset.fingerprint = data["fingerprint"]
    if "kmers" in data:
        kset.kmers = data["kmers"]
    if "counts" in data:
        kset.counts = data["counts"]
    return kset

def kmerSetFromHdf5(filePath):
    if not filePath.endswith(".hdf5"):
        filePath += ".hdf5"
    with h5py.File(filePath, 'r') as h5:
        assert h5.attrs["type"] == "KmerSet", "Not a KmerSet file!"
        kset = KmerSet(h5.attrs['k'])
        if "fingerprint" in h5:
            kset.fingerprint = np.array(h5["fingerprint"])
        if "kmers" in h5:
            kset.kmers = np.array(h5["kmers"])
        if "counts" in h5:
            kset.counts = np.array(h5["counts"])
    return kset



def kmerSetFromFile(filePath, k = DEFAULT_K):
    if filePath.endswith(".npz"):
        return kmerSetFromNpz(filePath, k)
    else:
        return kmerSetFromHdf5(filePath)


def similarityScore(kmers1, kmers2, scoring="jaccard"):
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
        score = intersection / ((kmers1.size + kmers2.size) / 2.0)
    elif scoring == "maxsize":
        # Use intersection / max_size (proper subset scores min/max)
        score = intersection / max(kmers1.size, kmers2.size)
    elif scoring == "reference":
        # Use intersection / size of reference (useful for comparing reads to assembled references)
        score = intersection / kmers2.size
    else:
        assert scoring in ("jaccard", "minsize", "maxsize", "meansize", "reference"), "unknown scoring method"
    return score


def similarityNumeratorDenominator(kmers1, kmers2, scoring="jaccard"):
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
        # Use intersection / size of reference (useful for comparing reads to assembled references)
        denom = kmers2.size
    else:
        assert scoring in ("jaccard", "minsize", "maxsize"), "unknown scoring method"
    return intersection, denom


class KmerSet:
    """
    Holds array of kmers and their associated counts & stats.
    """

    def __init__(self, k = DEFAULT_K):
        self.k = k
        # data arrays
        self.kmers = None
        self.counts = None
        self.fingerprint = None
        # stats
        self.nSeqs = 0
        self.nBases = 0
        self.nKmers = 0

    def __eq__(self, other):
        return self.k == other.k and np.array_equal(self.fingerprint, other.fingerprint) \
               and np.array_equal(self.kmers, other.kmers) and np.array_equal(self.counts, other.counts)

    def kmerizeFile(self, fileName, batchSize = 100000000, verbose = True, prune=0):
        seqFile = openSeqFile(fileName)
        batch = np.empty(batchSize, dtype=np.uint64)

        nSeqs = 0
        nBases = 0
        nKmers = 0
        nBatch = 0 # kmers in this batch
        pruned = False

        for seq in seqFile:
            nSeqs += 1
            seqLength = len(seq.seq)
            nBases += seqLength
            if nKmers + seqLength > batchSize:
                self.processBatch(batch, nSeqs, nBases, nKmers, verbose)
                if prune and self.singletons > prune:
                    self.pruneSingletons(verbose)
                    pruned = True
                nSeqs = 0
                nBases = 0
                nKmers = 0
            nKmers += kmerizer.kmerize_into_array(self.k, str(seq.seq), batch, nKmers)
        seqFile.close()
        self.processBatch(batch, nSeqs, nBases, nKmers, verbose)
        if pruned:
            self.pruneSingletons(verbose)

    def kmerizeSeq(self, seq):
        kmers = kmerizer.kmerize(self.k, seq)
        self.nSeqs += 1
        self.nBases += len(seq)
        self.nKmers = kmers.size
        self.kmers, self.counts = np.unique(kmers, return_counts=True)

    def processBatch(self, batch, nseqs, nbases, nkmers, verbose):
        self.nSeqs += nseqs
        self.nBases += nbases
        self.nKmers += nkmers

        newKmers, newCounts = np.unique(batch[:nkmers], return_counts=True)

        if self.kmers is None:
            self.kmers = newKmers
            self.counts = newCounts
        else:
            self.kmers, self.counts = kmerizer.merge_counts(self.kmers, self.counts, newKmers, newCounts)

        self.singletons = np.count_nonzero(self.counts == 1)
        if verbose:
            self.printStats()

    def pruneSingletons(self, verbose=False):
        keepers = self.counts > 1
        self.kmers = self.kmers[keepers]
        self.counts = self.counts[keepers]
        if verbose:
            print 'Pruned singletons:', self.kmers.size, 'distinct kmers remain'

    def mergeKmerSet(self, other):
        """Create new KmerSet by merging this with another"""
        newSet = KmerSet(self.k)
        newSet.kmers, newSet.counts = kmerizer.merge_counts(self.kmers, self.counts, other.kmers, other.counts)
        return newSet

    def printStats(self):
        print 'Seqs:', self.nSeqs, 'Bases:', self.nBases, 'Kmers:', self.nKmers, \
            'Distinct:', self.kmers.size, 'Singletons:', self.singletons

    def minHash(self, frac = 0.002):
        nkmers = int(round(self.kmers.size * frac))
        order = kmerizer.fnvhash_kmers(self.k, self.kmers).argsort()[:nkmers]
        self.fingerprint = self.kmers[order]
        self.fingerprint.sort()
        return self.fingerprint

    def freqFilter(self, minFreq = 1, maxFreq = None):
        condition = (self.counts >= minFreq)
        if maxFreq:
            condition &= (self.counts <= maxFreq)
        self.kmers = self.kmers[condition]
        self.counts = self.counts[condition]

    def spectrum(self):
        return np.unique(self.counts, return_counts=True)

    def spectrumMinMax(self, delta = .5, maxCopyNumber = 20):
        freq, counts = self.spectrum()
        minIndex = 0
        maxIndex = 0
        haveMin = False
        haveMax = False
        lastFreq = 0
        for i in xrange(freq.size):
            count = counts[i]
            zero = freq[i] > lastFreq + 1
            if haveMax and (zero or freq[i] > freq[maxIndex] * maxCopyNumber):
                break
            if haveMin:
                if count > counts[maxIndex]:
                    maxIndex = i
                if count < counts[maxIndex] * (1 - delta):
                    haveMax = True
            elif count > 1000 and count > counts[minIndex] * (1 + delta):
                haveMin = True
            elif zero or count < counts[minIndex]:
                minIndex = i
                maxIndex = i
            elif counts[i] < counts[minIndex]:
                minIndex = maxIndex = i
            lastFreq = freq[i]
            #print i, freq[i], zero, count, haveMin, freq[minIndex], haveMax, freq[maxIndex]
        if minIndex and maxIndex and counts[maxIndex] > counts[minIndex] * (1 + delta):
            return (freq[minIndex], freq[maxIndex], freq[i-1])
        return None

    def spectrumFilter(self, maxCopyNumber = 20):
        thresholds = self.spectrumMinMax()
        if thresholds:
            self.freqFilter(thresholds[0], thresholds[2])
        return thresholds


    def plotSpectrum(self, fileName = None, maxFreq = None):
        # to get kmer profile, count the counts!
        spectrum = self.spectrum()
        plt.semilogy(spectrum[0], spectrum[1])
        plt.grid = True
        if maxFreq:
            plt.xlim(0, maxFreq)
        plt.xlabel("Kmer Frequency")
        plt.ylabel("Number of Kmers")
        if fileName:
            plt.savefig(fileName)
        else:
            plt.show()

    def writeHistogram(self, fileName):
        spectrum = self.spectrum()
        with open(fileName, 'w') as hist:
            for i in xrange(spectrum[0].size):
                print >> hist, "%d\t%d" % (spectrum[0][i], spectrum[1][i])

    def save_npz(self, fileName, compress = False):
        kwargs = {'kmers': self.kmers, 'counts': self.counts}
        if self.fingerprint is not None:
            kwargs['fingerprint'] = self.fingerprint
        if compress:
            func = np.savez_compressed
        else:
            func = np.savez
        func(fileName, **kwargs)


    def save(self, fileName, compress = None, npz = False):
        """Save in HDF5 file format"""
        if npz:
            self.save_npz(fileName, compress)
            return
        if compress is True:
            compress = "gzip"
        if not fileName.endswith(".hdf5"):
            fileName += ".hdf5"
        with h5py.File(fileName, 'w') as h5:
            h5.attrs["type"] = np.string_("KmerSet")
            h5.attrs["k"] = self.k;
            if self.fingerprint is not None:
                h5.create_dataset("fingerprint", data=self.fingerprint, compression=compress)
            if self.kmers is not None:
                h5.create_dataset("kmers", data=self.kmers, compression=compress)
            if self.counts is not None:
                h5.create_dataset("counts", data=self.counts, compression=compress)

    def load_npz(self, fileName):
        npData = np.load(fileName)
        if 'kmers' in npData.files:
            self.kmers = npData['kmers']
        if 'counts' in npData.files:
            self.counts = npData['counts']
            self.nKmers = self.counts.sum()
        if 'fingerprint' in npData.files:
            self.fingerprint = npData['fingerprint']

    def load(self, fileName):
        if fileName.endswith(".npz"):
            self.load_npz(fileName)
            return
        with h5py.File(fileName, 'r') as h5:
            assert h5.attrs["type"] == "KmerSet", "Not a KmerSet file!"
            self.k = KmerSet(h5.attrs['k'])
            if "fingerprint" in h5:
                self.fingerprint = np.array(h5["fingerprint"])
            if "kmers" in h5:
                self.kmers = np.array(h5["kmers"])
            if "counts" in h5:
                self.counts = np.array(h5["counts"])










