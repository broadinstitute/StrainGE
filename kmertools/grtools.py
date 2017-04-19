#!/usr/bin/env python
"""Shared code for genome recovery tools"""

import os
import gzip
import cPickle
import numpy as np

min_qual = None
min_mq = None
min_confirm = None
min_gap = None
consensus = None
verbose = False

bases = "ACGT"

class Pileup:
    """
    Class to process pileup information; that is,
    all the alignment information corresponding
    to a given reference coordinate.
    """

    def __init__(self, chrom, scaffold, pos, pileup_reads = None):
        """
        :param scaffold: reference sequence
        :param pos: coordinate in reference (0-based)
        :param pileup_reads: pileup objects from pysam
        """
        self.chrom = chrom
        self.pos = pos
        self.reads = {}
        self.count = 0
        self.unmapped = 0
        self.bad = 0
        self.mq_total = 0
        self.mq_ss = 0
        self.qual_total = 0
        self.qual_ss = 0
        self.ref_count = 0
        self.ref_qual = 0
        self.ref_qual_ss = 0
        self.ref_mq = 0
        self.ref_mq_ss = 0
        self.hc = False
        # reads with other than reference base; list of tuples (base, qual, mq)
        self.others = {}

        self.refbase = scaffold[pos]
        if not pileup_reads:
            print "No reads for pileup"
            return
        for read in pileup_reads:
            self.add_read(read)

    def add_read(self, read):
        """Add a pysam.PileupRead object to this pileup"""
        alignment = read.alignment

        # if this is a paired read, make sure the pairs are properly aligned
        if (not alignment.is_paired) or (not alignment.is_proper_pair):
            self.bad += 1
            return

        # restrict ourselves to full-length alignments (not clipped)
        if alignment.query_alignment_length != alignment.query_length:
            # alignment is clipped
            self.bad += 1
            return

        # check that inferred insert size is at least read length
        tlen = alignment.template_length
        if abs(tlen) < alignment.query_length:
            self.bad += 1
            return

        # get base quality (note this is next base if deletion)
        pos = read.query_position_or_next
        qual = alignment.query_qualities[pos]
        if qual < min_qual:
            self.bad += 1
            return

        # base call must be real base (e.g., not N)
        base = alignment.query_sequence[pos]
        index = bases.find(base)
        if index < 0:
            self.bad += 1
            return

        # check for decent mapping quality
        mq = alignment.mapping_quality
        if mq < min_mq:
            self.unmapped += 1
            return

        # We're good! Update the pileup stats...
        self.count += 1
        self.mq_ss += mq**2
        self.mq_total += mq
        self.qual_ss += qual**2
        self.qual_total += qual
        # keep track of the reads in this pileup and the position in that read
        

        # read_name = alignment.query_name
        # if alignment.is_read1:
        #     read_name += ".1"
        # else:
        #     read_name += ".2"
        
        # self.reads[read_name] = pos
        if read.is_del:
            # using N as marker for deletion...
            # workaround until we can write actual deletions into vcf format
            self.add_other("N", qual, mq)
            #self.reads[read_name] = (pos, 'del')
            
        else:
            #self.reads[read_name] = (pos, base)
            if base == self.refbase:
                self.ref_count += 1
                self.ref_qual_ss += qual**2
                self.ref_qual += qual
                self.ref_mq_ss += mq**2
                self.ref_mq += mq
                
            else:
                self.add_other(base, qual, mq)
                if verbose:
                    print base, qual, mq

    def add_other(self, base, bq, mq):
        """
        Keep information about bases in pileup which differ from reference.
        :param base: alternate base or 'del' for deletion
        :param bq: base quality
        :param mq: mapping quality
        :return:
        """
        if not self.others:
            self.others = {}
        if base in self.others:
            (count, qss, mqss, qual, mapqual) = self.others[base]
        else:
            count = 0
            qual = 0
            mapqual = 0
            qss = 0
            mqss = 0
        count += 1
        qual += bq
        mapqual += mq
        qss += bq**2
        mqss += mq**2
        self.others[base] = (count, qss, mqss, qual, mapqual)

    def covered(self):
        """Does this pileup have enough data to consider this locus covered?"""
        # this seems wrong
        #return self.ref_qual >= min_confirm and self.count > self.bad
        
        # this seems better...
        #return self.qual_total >= min_confirm and self.count > (self.bad + self.unmapped)

        # why do we care if there are also bad reads as long as there are good ones???
        return self.qual_total >= min_confirm

    def unmappable(self):
        """Good quality sequence, but can't be mapped"""
        return self.unmapped > (self.bad + self.count)

    def ref_fraction(self):
        """Fraction of evidence which supports reference base"""
        if self.qual_total:
            return float(self.ref_qual) / self.qual_total
        else:
            return 0

    def confirmed(self):
        """Does this pileup confirm the reference?"""
        return self.ref_qual >= min_confirm and \
               (self.ref_qual == self.qual_total or self.ref_fraction() > consensus)
    
    def _get_best_snp(self):
        """Returns the SNP with the highest sum of RMS base quality and mapping quality, if two are qual, return higher count. All else equal, alphabetical."""
        best_snp = None
        best_score = 0
        best_count = 0
        if len(self.others) == 1:
            return self.others.keys()[0]
        for snp in self.others:
            # score is just sum of RMS(BQ) and RMS(MQ)
            score = np.sqrt(float(self.others[snp][1])/self.others[snp][0]) + np.sqrt(float(self.others[snp][2])/self.others[snp][0])
            count = self.others[snp][0]
            if score > best_score:
                best_snp = snp
                best_score = score
                best_count = count
            elif score == best_score:
                if count > best_count:
                    best_snp = snp
                    best_score = score
                    best_count = count
                elif count == best_count:
                    if verbose: print "Warning: multiple SNPs have same score and count!"
                    return False
        
        return best_snp
        
        #return max(self.others.keys(), key = lambda snp: (np.sqrt(float(self.others[snp][1])/self.others[snp][0])+np.sqrt(float(self.others[snp][2])/self.others[snp][0]), self.others[snp][0], snp))
    
    def sort_alts(self):
        return sorted(self.others.keys(), key = lambda snp: (-np.sqrt(float(self.others[snp][1])/self.others[snp][0]) + np.sqrt(float(self.others[snp][2])/self.others[snp][0]), -self.others[snp][0], snp))
    
    def base_call(self):
        """Call another base...just print out stats for now"""
        rf = self.ref_fraction()
        if rf < consensus:
            if verbose: print 'SNP?', self.pos, self.refbase, self.count, self.bad, rf, self.others
            
            if len(self.others) == 0:
                if verbose: print "No evidence of SNPs"
            else:
                # get SNP with highest proportion if >1 SNP
                base = self._get_best_snp()
                #if base and self.others[base][0] >= (consensus * self.count) and self.others[base][3] >= min_confirm:
                if base and (float(self.others[base][3]) / self.qual_total) >= consensus and self.others[base][3] >= min_confirm:
                    # 90% of base quality matches SNP
                    if verbose: print "SNP confirmed %s" % base
                    return base
                else:
                    # no consensus
                    if verbose: print "No confirmed SNP"
        
        return None
    
    def high_coverage(self, high):
        """Flag pileup for too much coverage"""
        self.hc = self.count > high
    
    def __str__(self):
        return "<Pileup n=%d rc=%d bad=%d rq=%d/%d mq=%d/%d o=%s>" % (self.count, self.ref_count, self.bad,
                                                                      self.ref_qual, self.qual_total,
                                                                      self.ref_mq, self.mq_total, str(self.others))

class Pileups:
    def __init__(self, pileups=None, name=None, keep_covered=True):
        self.pileups = {}
        self.reads = {}
        self.name = name
        self.length = 0

        for pileup in pileups:
            if verbose:
                print "Adding pileup", str(pileup)
            self.add_pileup(pileup, keep_covered=keep_covered)
            self.length += 1

    def add_pileup(self, pileup, keep_covered=False):
        """Add pileups based on scaffold and position"""

        if keep_covered and not pileup.covered():
            return

        chrom = pileup.chrom
        if chrom not in self.pileups:
            self.pileups[chrom] = {}
        pos = pileup.pos
        self.pileups[chrom][pos] = pileup

        # store reads and their positions, match back to pileup
        for read in pileup.reads:
            if read not in self.reads:
                self.reads[read] = {}
            if pileup.reads[read][1] == 'del':
                continue
            readpos = pileup.reads[read][0]
            if readpos in self.reads[read]:
                print "Duplicate read position found", read, readpos
                print self.reads[read][readpos]
                print chrom, pos
            self.reads[read][readpos] = (chrom, pos)

    def __len__(self):
        return self.length

    def __str__(self):
        return "<Pileups Class> containing %i pileups" % self.length


def save_pileups(pileups, out):
    """Save all pileups to a pickled file"""

    with gzip.open(out, 'wb') as w:
        cPickle.dump(pileups, w)

def load_pileups(pkl_file, keep_covered=False):
    """Load straingr pileups"""
    with gzip.open(pkl_file) as f:
        pileups = Pileups(cPickle.load(f), name=os.path.basename(pkl_file), keep_covered=keep_covered)

        return pileups
