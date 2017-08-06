#!/usr/bin/env python

import sys
import re
import kmertools
import shelve
import numpy as np

class BarcodeIndex:
    """
    Holds index mapping read names to barcodes. Index stored in array of two columns,
    where first holds hash of read name and second has binary encoded barcode.
    """

    BLOCK_SIZE = 1000000
    BASES = "ACGT"
    BARCODE_LENGTH = 16

    def __init__(self, filename=None):
        self.count = 0
        self.size = 0
        self.bcindex = None
        self.sorted = False
        if filename:
            self.load(filename)
            self.sorted = True

    def add(self, readname, barcode):
        if self.count >= self.size:
            self.extend()
            print(self.count)
        self.bcindex[self.count, 0] = hash(readname)
        self.bcindex[self.count, 1] = self.barcode_to_int64(barcode)
        self.count += 1
        self.sorted = False

    def extend(self):
        """
        Extend index array by BLOCK_SIZE.
        """
        block = np.zeros((self.BLOCK_SIZE, 2), dtype=np.int64)
        self.bcindex = np.concatenate((self.bcindex, block)) if self.bcindex is not None else block
        self.size = len(self.bcindex)

    def barcode_to_int64(self, barcode):
        """
        :param barcode: 10X barcode bases
        :return: binary encoding using 2 bits per base (A=0,C=1,G=2,T=3)
        """
        value = 0
        for c in barcode:
            n = self.BASES.index(c)
            value = (value << 2) | n
        return value

    def int64_to_barcode(self, value):
        """
        :param value: binary
        :return: string of barcode bases
        """
        bases = [self.BASES[(value >> (2 * (self.BARCODE_LENGTH - 1 - i))) & 3] for i in xrange(self.BARCODE_LENGTH)]
        return "".join(bases)

    def lookup(self, readname):
        """
        Find barcode of given read name
        :param readname: read name string
        :return: barcode string
        """
        if not self.sorted:
            self.sort()
        h = hash(readname)
        n = np.apply_along_axis(lambda a: a.searchsorted(h), axis=0, arr=self.bcindex)[0]
        print(n, self.bcindex)
        if self.bcindex[n, 0] == h:
            return self.int64_to_barcode(self.bcindex[n, 1])
        return None

    def sort(self):
        self.bcindex = self.bcindex[:self.count]
        self.bcindex = self.bcindex[self.bcindex[:, 0].argsort()]
        self.sorted = True

    def save(self, filename):
        """
        Save barcode index array to file
        :param filename: Filename for compressed npz file to hold 'bcindex' array
        """
        if not self.sorted:
            self.sort()
        np.savez_compressed(filename, bcindex=self.bcindex)

    def load(self, filename):
        """
        Restore barcode index array from file
        :param filename: Filename of compressed npz file containing 'bcindex' array
        """
        data = np.load(filename)
        self.bcindex = data["bcindex"]
        print(self.bcindex)
        self.count = len(self.bcindex)
        self.size = self.count



tagre = re.compile("BX:Z:([ACGT]+)-[0-9]+")
fastq = kmertools.openSeqFile(sys.argv[1])

barcode_index = BarcodeIndex()

for fq in fastq:
    match = tagre.search(fq.description)
    if match:
        barcode = match.groups()[0]
    else:
        continue
    barcode_index.add(fq.name, barcode)

barcode_index.save("bcindex.npz")
