import sys
import gzip
import bz2
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import kmerizer

A = 0
C = 1
G = 2
T = 3

BASES = "ACGT"

# A random 64-bit number used in hashing function
HASH_BITS = 0x29679e096c8c07bf

def openSeqFile(fileName):
    """
    Open a sequence file with SeqIO; can be fasta or fastq with optional gz or bz2 compression.
    Assumes fasta unless ".fastq" or ".fq" in the file name.
    :param fileName:
    :return: SeqIO.parse object
    """

    components = fileName.split('.')
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

def loadKmers(fileName):
    return loadNpz(fileName, "kmers")

def loadCounts(fileName):
    return loadNpz(fileName, "counts")

def loadFingerprint(fileName):
    return loadNpz(fileName, "fingerprint")

class KmerSet:
    """
    Holds array of kmers and their associated counts & stats.
    """

    def __init__(self, k = 31):
        self.k = k
        # data arrays
        self.kmers = None
        self.counts = None
        self.fingerprint = None
        # stats
        self.nSeqs = 0
        self.nBases = 0
        self.nKmers = 0

    def kmerizeFile(self, fileName, batchSize = 100000000, verbose = True):
        seqFile = openSeqFile(fileName)
        batch = np.empty(batchSize, dtype=np.uint64)

        nSeqs = 0
        nBases = 0
        nKmers = 0
        nBatch = 0 # kmers in this batch

        for seq in seqFile:
            nSeqs += 1
            seqLength = len(seq)
            nBases += seqLength
            if nKmers + seqLength > batchSize:
                self.processBatch(batch, nSeqs, nBases, nKmers, verbose)
                nSeqs = 0
                nBases = 0
                nKmers = 0
            nKmers += kmerizer.kmerize_into_array(self.k, str(seq.seq), batch, nKmers)
        seqFile.close()
        self.processBatch(batch, nSeqs, nBases, nKmers, verbose)

    def processBatch(self, batch, nseqs, nbases, nkmers, verbose):
        self.nSeqs += nseqs
        self.nBases += nbases
        self.nKmers += nkmers

        newKmers, newCounts = np.unique(batch[:nkmers], return_counts=True)

        if type(self.kmers) == type(None):
            self.kmers = newKmers
            self.counts = newCounts
        else:
            self.kmers, self.counts = kmerizer.merge_counts(self.kmers, self.counts, newKmers, newCounts)

        if verbose:
            self.printStats()

    def printStats(self):
        print 'Seqs:', self.nSeqs, 'Bases:', self.nBases, 'Kmers:', self.nKmers, \
            'Distinct:', self.kmers.size, 'Singletons:', np.count_nonzero(self.counts == 1)

    def hashKmers(self, kmers):
        mask = (1 << (2 * self.k)) - 1
        hashedKmers = kmers ^ (HASH_BITS & mask)
        hashedKmers.sort()
        return hashedKmers

    def unHashKmers(self, kmers):
        # un-do the hash function; this trival one self-reverses
        return self.hashKmers(kmers)

    def minHash(self, nkmers = 10000):
        self.fingerprint = self.hashKmers(self.kmers)[:nkmers]
        self.fingerprint = self.unHashKmers(self.fingerprint)
        return self.fingerprint

    def freqFilter(self, minFreq = 1, maxFreq = None):
        condition = (self.counts >= minFreq)
        if maxFreq:
            condition &= (self.counts <= maxFreq)
        self.kmers = self.kmers[condition]
        self.counts = self.counts[condition]

    def spectrum(self):
        return np.unique(self.counts, return_counts=True)

    def spectrumMinMax(self, delta = .5):
        freq, counts = self.spectrum()
        minIndex = 0
        maxIndex = 0
        for i in xrange(freq.size):
            count = counts[i]
            if minIndex and counts[i] > counts[maxIndex]:
                maxIndex = i
            elif counts[i] < counts[minIndex]:
                minIndex = maxIndex = i
            if minIndex and maxIndex and counts[i] < counts[maxIndex] * (1 - delta):
                return (freq[minIndex], freq[maxIndex])
        return None

    def spectrumFilter(self, maxCopyNumber = 20):
        thresholds = self.spectrumMinMax()
        if thresholds:
            self.freqFilter(thresholds[0], thresholds[1] * maxCopyNumber)
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

    def save(self, fileName, compress = False):
        kwargs = {'kmers': self.kmers, 'counts': self.counts}
        if type(self.fingerprint) != type(None):
            kwargs['fingerprint'] = self.fingerprint
        if compress:
            func = np.savez_compressed
        else:
            func = np.savez
        func(fileName, **kwargs)

    def load(self, fileName):
        npData = np.load(fileName)
        if 'kmers' in npData.files:
            self.kmers = npData['kmers']
        if 'counts' in npData.files:
            self.counts = npData['counts']
            self.nKmers = self.counts.sum()
        if 'fingerprint' in npData.files:
            self.fingerprint = npData['fingerprint']

    def kmerString(self, kmer):
        return ''.join([BASES[(kmer >> shift) & 3] for shift in xrange(2 * self.k - 2, -1, -2)])










