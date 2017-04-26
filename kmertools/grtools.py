#!/usr/bin/env python
"""Shared code for genome recovery tools"""

import os
import math
import pysam
import gzip
import cPickle
import numpy as np

min_qual = 5
min_mq = 5
min_confirm = 50
consensus = 0.9
min_gap = 2000
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
        
        read_name = alignment.query_name
        if alignment.is_read1:
            read_name += ".1"
        else:
            read_name += ".2"
        
        # self.reads[read_name] = pos
        if read.is_del:
            # using N as marker for deletion...
            # workaround until we can write actual deletions into vcf format
            self.add_other("N", qual, mq)
            self.reads[read_name] = (pos, 'del')
            
        else:
            self.reads[read_name] = (pos, base)
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
                snp_base = self._get_best_snp()
                #if base and self.others[base][0] >= (consensus * self.count) and self.others[base][3] >= min_confirm:
                if snp_base and (float(self.others[snp_base][3]) / self.qual_total) >= consensus and self.others[snp_base][3] >= min_confirm:
                    # 90% of base quality matches SNP
                    if verbose: print "SNP confirmed %s" % snp_base
                    return snp_base
                else:
                    # no consensus
                    if verbose: print "No confirmed SNP"
        
        return None
    
    def high_coverage(self, high):
        """Flag pileup for too much coverage"""
        self.hc = self.count > high
        if self.hc:
            return 1
        else:
            return 0
    
    def __str__(self):
        return "<Pileup n=%d rc=%d bad=%d rq=%d/%d mq=%d/%d o=%s>" % (self.count, self.ref_count, self.bad,
                                                                      self.ref_qual, self.qual_total,
                                                                      self.ref_mq, self.mq_total, str(self.others))

class Pileups:
    """Class of Pileups and reads for straingr"""
    def __init__(self, reference, bamfile):
        self.pileups = {}
        self.reads = {}
        self.nPileups = 0
        self.length = 0
        self.confirmed = 0
        self.covered = 0
        self.snps = 0
        self.unmapped = 0
        self.goodcoverage = 0
        self.highcoverage = 0

        print "Scanning BAM file: %s" % bamfile
        bam = pysam.AlignmentFile(bamfile, "rb")
    
        for scaffold in bam.references:
            try:
                self.process_scaffold(bam, scaffold, reference[scaffold])
            except KeyboardInterrupt:
                raise e
            except Exception as e:
                print "Error in straingr: %s" % e
                continue
    
    def __len__(self):
        return self.nPileups

    def process_scaffold(self, bam, scaffold, refseq):
        """Scan the pileups for each locus in the scaffold"""

        self.pileups[scaffold] = {}
        length = len(refseq)
        print "Processing", scaffold, length
        confirmed = 0
        covered = 0
        snps = 0
        unmapped = 0
        goodcoverage = 0
        highcoverage = 0

        last_covered = -1
        gaps = []
        
        for column in bam.pileup(scaffold):
            self.nPileups += 1
            refpos = column.reference_pos
            pileup = Pileup(column.reference_name, refseq, refpos, column.pileups)

            for read in pileup.reads:
                if read not in self.reads:
                    self.reads[read] = {}
                (pos, base) = pileup.reads[read]
                self.reads[read][pos] = (base, scaffold, refpos)
            
            # don't keep this info, as it's redundant
            del pileup.reads

            if verbose:
                refbase = refseq[refpos]
                print "Ref:", column.reference_name, refpos, refbase, column.nsegments
            goodcoverage += pileup.count
            if pileup.covered():
                covered += 1
                if pileup.confirmed():
                    confirmed += 1
                elif pileup.base_call():
                    snps += 1
                if refpos - last_covered > min_gap:
                    gap = (last_covered + 1, refpos - last_covered)
                    print "Coverage gap:", gap[0], gap[1]
                    gaps.append(gap)
                last_covered = refpos
            elif pileup.unmappable():
                unmapped += 1
                # not a real gap, just can't map to this region
                if refpos - last_covered > min_gap:
                    gap = (last_covered + 1, refpos - last_covered)
                    print "Coverage gap:", gap[0], gap[1]
                    gaps.append(gap)
                last_covered = refpos
            
            if verbose:
                print pileup, pileup.confirmed()
            self.pileups[scaffold][refpos] = pileup

        coverage = float(goodcoverage) / float(length)
        mixed = covered - (confirmed + snps)

        print "good coverage: %.1fx" % (coverage,)
        print "covered: %d %.1f%%" % (covered, pct(covered, length))
        print "confirmed: %d %.2f%%" % (confirmed, pct(confirmed, covered))
        print "snps: %d %.3f%%" % (snps, pct(snps, covered))
        if snps > 0:
            print "snp rate: %.0f" % (float(covered) / float(snps))
        print "mixed: %d %.3f%%" % (mixed, pct(mixed, covered))
        if mixed > 0:
            mixed_rate = float(covered) / float(mixed)
            if mixed_rate > 0:
                mixed_quality = math.log10(mixed_rate) * 10.0
            else:
                mixed_quality = 0
            print "mixed rate: %.0f Q%.0f" % (mixed_rate, mixed_quality)
        print "gaps:", len(gaps), "totaling", sum([g[1] for g in gaps])
        print "unmapped: %d %.1f%%" % (unmapped, pct(unmapped, length))
        
        # keep track of total values in the Pileups class variables
        self.length += length
        self.confirmed += confirmed
        self.covered += covered
        self.snps += snps
        self.unmapped += unmapped
        self.goodcoverage += goodcoverage

        if len(self.pileups[scaffold]):
            avg_count = goodcoverage / len(self.pileups[scaffold])
            if avg_count > 3:
                # threshold is 99.9999%ile of poissons distribution at average coverage
                threshold = int(np.percentile(np.random.poisson(avg_count, self.nPileups), 99.9999))
            else:
                # at lower coverages, just set it to 15
                threshold = 15
            for refpos in self.pileups[scaffold]:
                highcoverage += self.pileups[scaffold][refpos].high_coverage(threshold)
            
            print "Abnormally high coverage: %d %.1f%% (expect 0.01%% false positive)" % (highcoverage, pct(highcoverage, length))

            self.highcoverage += highcoverage

        

        


def pct(numerator, denominator):
    """Makes into a percent"""
    if numerator > 0 and denominator > 0:
        return (100.0 * numerator) / denominator
    return 0.0


def save_pileups(pileups, out):
    """Save all pileups to a pickled file"""

    with gzip.open(out, 'wb') as w:
        cPickle.dump(pileups, w)

def load_pileups(pkl_file):
    """Load saved pileups from pickled file"""

    with gzip.open(pkl_file, 'rb') as f:
        pileups = cPickle.load(f)
        # if type(pileups) is not Pileups:
        #     print 'Not a valid pickled pileups file'
        #     return
        return pileups

# def load_pileups(pkl_file, keep_covered=False):
#     """Load straingr pileups"""
#     with gzip.open(pkl_file) as f:
#         pileups = Pileups(cPickle.load(f), name=os.path.basename(pkl_file), keep_covered=keep_covered)

#         return pileups
