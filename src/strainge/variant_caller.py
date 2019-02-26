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

import math
import logging
import functools
from enum import IntFlag, auto
from typing import Dict  # noqa

import numpy
from scipy.stats import poisson

from strainge import kmertools, utils
from strainge.utils import pct

logger = logging.getLogger(__name__)


class Allele(IntFlag):
    """Enum for possible alleles at a position. Derives from `enum.IntFlag` so
    values can be combined like this:

    >>> Allele.A | Allele.T
    <Allele.A|T: 9>

    This is useful to indicate that multiple alleles are present at a given
    genomic location.
    """

    N = 0
    A = auto()
    C = auto()
    G = auto()
    T = auto()
    INS = auto()
    DEL = auto()

    @classmethod
    def from_str(cls, base):
        """Create a new `Allele` object from a single character string."""

        if base not in cls.__members__:
            return Allele.N

        return cls.__members__[base]

    def __iter__(self):
        for allele in Allele:
            if self.value & allele:
                yield allele

    def __str__(self):
        alleles = list(self)
        if len(alleles) == 1:
            return _allele_to_str(alleles[0])
        else:
            return ",".join(str(v) for v in self)


@functools.lru_cache(maxsize=8)
def _allele_to_str(value):
    rev_mapping = {v: k for k, v in Allele.__members__.items()}

    return rev_mapping[value]


ALLELE_MASKS = numpy.array([v for v in Allele if v != Allele.N])

ALLELE_INDEX = {
    a: i for i, a in enumerate(v for v in Allele if v != Allele.N)
}


def poisson_coverage_cutoff(mean, cutoff=0.9999999):
    """
    Calculate the Poisson CDF and find where it reaches the cutoff. For
    higher coverages, use linear instead.

    Default cutoff is one part in 10M, so not likely to occur in a typical
    bacterial genome.
    """

    if mean < 50:
        # for lower coverage, use poisson
        return poisson.ppf(cutoff, mean)
    else:
        return int(math.ceil(mean * 1.5 + 15.0))


def scale_min_gap_size(min_gap, mean_coverage):
    """
    Attempt to scale the minimum significant uncovered region by coverage
    (for low coverage)

    :param min_gap: original min gap
    :param mean_coverage: expected coverage
    :return: scaled minimum gap
    """
    lw = utils.lander_waterman(mean_coverage)

    return int(min_gap / lw) if lw > 0 else min_gap


class Reference:
    """
    Holds reference sequence information...scaffolds in Biopython SeqRecord
    objects.
    """
    def __init__(self, fasta):
        self.fasta = fasta
        self.scaffolds = {
            scaffold.name: scaffold
            for scaffold in kmertools.open_seq_file(fasta)
        }
        self.lengths = [len(s) for s in self.scaffolds.values()]
        self.length = sum(self.lengths)

        logger.info("Reference %s has %d scaffolds with a total of %d bases.",
                    fasta, len(self.scaffolds), self.length)

    def scaffold_coord(self, coord):
        """
        Turn a zero-based genome-wide coordinate into a scaffold & coordinate
        within scaffold (1-based)

        :param coord: zero-based genome-wide coordinate
        :return: (scaffold, scaffoldCoord)
        """
        offset = 0
        for scaffold, length in zip(self.scaffolds.values(), self.lengths):
            if coord < offset + length:
                return scaffold.name, coord + 1 - offset
            offset += length

    def scaffold_to_genome_coord(self, scaffold_name, coord):
        """
        Turn a 1-based scaffold coordinate into a 0-based genome-wide
        coordinate.

        :param scaffold_name: scaffold name
        :param coord: 1-based scaffold coordinate
        :return: genomeCoord
        """
        offset = 0
        for scaffold, length in zip(self.scaffolds, self.lengths):
            if scaffold.name == scaffold_name:
                return offset + coord - 1
            offset += length

    def get_sequence(self, name, coord, length=1):
        return self.scaffolds[name].seq[coord-1:coord+length-1]


class VariantCallData:
    """
    This class holds all data and statistics needed for variant calling. The
    data is stored per contig/scaffold in the reference.
    """

    def __init__(self, reference, min_gap_size):
        """
        :param reference: The fasta file used as reference for read alignment.
        :type reference: Reference
        """

        self.reference = reference
        self.min_gap_size = min_gap_size

        self.scaffolds_data = {
            name: ScaffoldCallData(name, len(scaffold.seq))
            for name, scaffold in self.reference.scaffolds.items()
        }  # type: Dict[str, ScaffoldCallData]

        self.mean_coverage = 0.0
        self.median_coverage = 0

    def build_refmask(self):
        for name, scaffold in self.reference.scaffolds.items():
            logger.info("Building refmask for scaffold %s", name)

            for i, base in enumerate(scaffold.seq.upper()):
                scaffold.refmask[i] = Allele.from_str(base)

    def bad_read(self, scaffold, pos):
        self.scaffolds_data[scaffold].bad[pos] += 1

    def low_mapping_quality(self, scaffold, pos):
        if pos >= self.scaffolds_data[scaffold].length:
            logger.warning("Position %d for scaffold %s of length %d out of "
                           "bounds, ignoring!", pos, scaffold,
                           self.scaffolds_data[scaffold].length)
            return

        self.scaffolds_data[scaffold].lowmq_count[pos] += 1

    def update_mapping_quality(self, scaffold, pos, mapping_quality):
        self.scaffolds_data[scaffold].mq_sum[pos] += mapping_quality

    def good_read(self, scaffold, pos, allele, mapping_quality):
        ix = ALLELE_INDEX[allele]
        self.scaffolds_data[scaffold].alleles[pos, 0, ix] += 1
        self.scaffolds_data[scaffold].alleles[pos, 1, ix] += mapping_quality

    def analyze_coverage(self):
        for scaffold in self.scaffolds_data.values():
            scaffold.calculate_coverage()

        all_coverage = numpy.concatenate([s.coverage for s in
                                          self.scaffolds_data.values()])
        self.mean_coverage = numpy.sum(all_coverage) / self.reference.length
        self.median_coverage = numpy.median(all_coverage)

        return self

    def call_alleles(self, min_pileup_qual, min_qual_frac):
        for scaffold in self.scaffolds_data.values():
            scaffold.call_alleles(min_pileup_qual, min_qual_frac)

        return self

    def find_gaps(self):
        for scaffold in self.scaffolds_data.values():
            scaffold.find_gaps(self.min_gap_size)

        return self

    def summarize(self):
        """
        Summarize all earlier calculated statistics into a global overview for
        the whole genome.
        """

        total_callable = 0
        all_coverages = []
        total_confirmed = 0
        total_snps = 0
        total_multi = 0
        total_lowmq = 0
        total_high_cov = 0
        total_gaps = 0
        total_gap_length = 0

        for scaffold in self.scaffolds_data.values():
            # Locations with strong evidence for the reference base
            confirmed = (scaffold.strong & scaffold.refmask)

            # Locations where we have strong evidence something else than the
            # reference (not mutually exclusive with the above!)
            snps = (scaffold.strong & ~scaffold.refmask)

            # Locations where we have strong evidence for multiple bases (could
            # be both reference or not)
            # (x & (x - 1)) turns of highest bit in x, will result in a
            # non-zero value if we have multiple alleles (that is, multiple
            # bits set)
            multi = (scaffold.strong & (scaffold.strong - 1)) > 0

            # Consider a locus callable if we have a strong call
            num_callable = numpy.count_nonzero(scaffold.strong)
            callable_pct = pct(num_callable, scaffold.length)
            total_callable += num_callable

            # Don't take abnormal high coverage regions into account when
            # calculating mean coverage
            normal_coverage = ~scaffold.high_coverage
            summed_coverage = (scaffold.coverage[normal_coverage].sum() +
                               scaffold.lowmq.sum())
            coverage = summed_coverage / scaffold.length

            median_coverage = numpy.median(scaffold.coverage + scaffold.lowmq)
            all_coverages.append(coverage)

            num_confirmed = numpy.count_nonzero(confirmed)
            confirmed_pct = pct(num_confirmed, num_callable)
            total_confirmed += num_confirmed

            num_snps = numpy.count_nonzero(snps)
            snp_pct = pct(num_snps, num_callable)
            total_snps += num_snps

            num_multi = numpy.count_nonzero(multi)
            multi_pct = pct(num_multi, num_snps)
            total_multi += num_multi

            num_lowmq = numpy.count_nonzero(scaffold.lowmq)
            lowmq_pct = pct(num_lowmq, scaffold.length)
            total_lowmq += num_lowmq

            num_high_cov = numpy.count_nonzero(scaffold.high_coverage)
            high_pct = pct(num_high_cov, scaffold.length)
            total_high_cov += num_high_cov

            num_gaps = len(scaffold.gaps)
            gap_length = sum(g.length for g in scaffold.gaps)
            total_gaps += num_gaps
            total_gap_length += gap_length

            yield {
                "name": scaffold.name,
                "length": scaffold.length,
                "coverage": coverage,
                "median": median_coverage,
                "callable": num_callable,
                "callablePct": callable_pct,
                "confirmed": num_confirmed,
                "confirmedPct": confirmed_pct,
                "snps": num_snps,
                "snpPct": snp_pct,
                "multi": num_multi,
                "multiPct": multi_pct,
                "lowmq": num_lowmq,
                "lowmqPct": lowmq_pct,
                "high": num_high_cov,
                "highPct": high_pct,
                "gapCount": num_gaps,
                "gapLength": gap_length
            }

        # Return one last entry with all statistics for the genome as a whole
        yield {
            "name": "TOTAL",
            "length": self.reference.length,
            "coverage": (
                # Weigh by scaffold length
                sum(v * l for v, l in
                    zip(all_coverages, self.reference.lengths)) /
                self.reference.length
            ),
            "median": self.median_coverage,
            "callable": total_callable,
            "callablePct": pct(total_callable, self.reference.length),
            "confirmed": total_confirmed,
            "confirmedPct": pct(total_confirmed, total_callable),
            "snps": total_snps,
            "snpPct": pct(total_snps, total_callable),
            "multi": total_multi,
            "multiPct": pct(total_multi, total_snps),
            "lowmq": total_lowmq,
            "lowmqPct": pct(total_lowmq, self.reference.length),
            "high": total_high_cov,
            "highPct": pct(total_high_cov, self.reference.length),
            "gapCount": total_gaps,
            "gapLength": total_gap_length
        }


class ScaffoldCallData:
    """
    Contains statistics about pileups for a given scaffold.
    """

    def __init__(self, name, length):
        self.name = name
        self.length = length

        self.refmask = numpy.zeros((self.length,), dtype=numpy.uint8)

        # Store for each position and per possible allele the counts and sum
        # of base qualities. We use len(Allele)-1 because we store nothing
        # for Allele.N.
        self.alleles = numpy.zeros((self.length, 2, len(Allele)-1),
                                   dtype=numpy.uint32)

        # Number of reads rejected for some reason
        self.bad = numpy.zeros((self.length,), dtype=numpy.uint32)

        # Number of reads with low mapping quality (for example a repetitive
        # region)
        self.lowmq_count = numpy.zeros((self.length,), dtype=numpy.uint32)

        # Regions in this scaffold with more low mapping quality reads than
        # "good" reads
        self.lowmq = None

        # Sum of mapping qualities per position
        self.mq_sum = numpy.zeros((self.length,), dtype=numpy.uint32)

        # Alleles which we can call confidently, will be filled by
        # `call_alleles`.
        self.strong = None

        # Alleles with weak evidence
        self.weak = None

        # Coverage of good quality reads, will be computed from alleles and
        # low mapping quality count
        self.coverage = None

        # Positions marked as too high coverage (conserved genes for example)
        self.high_coverage = None

        self.mean_coverage = 0.0
        self.median_coverage = 0
        self.coverage_cutoff = 0

        self.gaps = []

    def calculate_coverage(self):
        """
        Calculate coverage for each position, which is calculated from
        observed good reads + the low mapping quality counts. Furthermore,
        we determine which positions have exceptionally high coverage.

        Based on coverage, set a high limit for plausible coverage of a
        given position based on Poisson probability; because we are
        dealing with metagenomic samples, too improbably high coverage is
        likely a result of reads from other organisms aligning to conserved
        regions, so we can't make a confident call on our target organism.
        Here, we use the median coverage rather than the mean since our
        coverage might be dominated by conserved regions.
        """
        self.coverage = self.alleles[:, 0].sum(axis=-1) + self.lowmq_count
        self.mean_coverage = numpy.sum(self.coverage) / self.length
        self.median_coverage = numpy.median(self.coverage)

        self.coverage_cutoff = poisson_coverage_cutoff(
            max(0.5, self.median_coverage))

        logger.info("Scaffold %s has mean coverage %.2f (median: %d). High "
                    "coverage cutoff: %d.", self.name, self.mean_coverage,
                    self.median_coverage, self.coverage_cutoff)

        self.high_coverage = self.coverage > self.coverage_cutoff

    def call_alleles(self, min_pileup_qual, min_qual_frac):
        quals = self.alleles[:, 1]
        qual_sums = quals.sum(axis=-1)
        qual_fraction = numpy.divide(quals, qual_sums[:, numpy.newaxis],
                                     where=qual_sums[:, numpy.newaxis] > 0)

        evidence = quals > 0
        # ALLELE_MASKS is an array with per allele its bit value.
        # By multiplying it with the above boolean array and summing the
        # result, we set each bit for each allele for which we have observed
        # evidence.
        self.weak = (evidence * ALLELE_MASKS[numpy.newaxis, :]).sum(axis=-1)

        confirmed = ((quals > min_pileup_qual) &
                     (qual_fraction > min_qual_frac))
        self.strong = (confirmed * ALLELE_MASKS[numpy.newaxis, :]).sum(axis=-1)

        # Remove any calls in too high coverage regions
        self.weak[self.high_coverage] = 0
        self.strong[self.high_coverage] = 0

    def find_gaps(self, min_size):
        """
        Find coverage gaps, indicating possible deleted portions in our
        strain vs the reference (e.g., possible recombination events). Gaps
        are defined here as insufficient evidence to call anything, taking
        into account regions of poor mapping quality won't be called as a
        matter of course (but may still be "covered").

        Only gaps larger than `min_size` are reported, and this value is
        scaled taking mean coverage of this scaffold into account.
        """
        min_size = scale_min_gap_size(min_size, self.mean_coverage)
        logger.info("%s: scaled min-gap size %.2f at mean coverage %.2f",
                    self.name, min_size, self.mean_coverage)

        depth = self.alleles[:, 0].sum(axis=-1)
        self.lowmq = ((self.lowmq_count > 1) & (self.lowmq_count > depth))

        # Covered is either: 1) we can make a strong call 2) we have low
        # mapping quality reads there (potentially from repetitive regions)
        covered_array = ((self.strong > 0) | self.lowmq)

        self.gaps = [
            group for group in utils.find_consecutive_groups(covered_array,
                                                             min_size)
            # group with all zeros is a region with no coverage
            if not numpy.any(group.data)
        ]

    def depth(self, loc):
        """
        :return: Count of all good reads
        :rtype: int
        """
        return self.alleles[loc, 0].sum()

    def qual_total(self, loc):
        """
        :return: Sum of all quality evidence
        :rtype: int
        """
        return self.alleles[loc, 1].sum()

    def total_depth(self, loc):
        """
        :return: Count of all reads, including those rejected for some reason
        :rtype: int
        """
        return self.depth(loc) + self.lowmq_count[loc]

    def ref_count(self, loc):
        return self.allele_count(loc, self.refmask[loc])

    def ref_qual(self, loc):
        """
        :return: sum of quality evidence for reference base (int)
        """
        ix = ALLELE_INDEX[self.refmask[loc]]
        return self.alleles[loc, 1, ix]

    def ref_fraction(self, loc):
        """
        :return: Fraction of evidence which supports reference base (float)
        """
        return self.ref_qual(loc) / self.qual_total(loc)

    def allele_count(self, loc, allele):
        return self.alleles[loc, 0, ALLELE_INDEX[allele]]

    def allele_qual(self, loc, allele):
        return self.alleles[loc, 1, ALLELE_INDEX[allele]]


class VariantCaller:
    """
    This class collects read alignments and updates any call statistics per
    scaffold.
    """

    def __init__(self, min_qual, min_pileup_qual, min_qual_frac,
                 min_mapping_quality, min_gap_size):
        self.min_qual = min_qual
        self.min_pileup_qual = min_pileup_qual
        self.min_qual_frac = min_qual_frac
        self.min_mapping_quality = min_mapping_quality
        self.min_gap_size = min_gap_size

    def process(self, reference, pileup_iter):
        """
        Process the pileups from a BAM file and collect all statistics and
        data reequired for variant calling
        :param reference: Which reference used for alignment
        :type reference: Reference
        :param pileup_iter: Iterable of pysam.PileupColumn objects
        :return:
        """
        call_data = VariantCallData(reference, self.min_gap_size)
        call_data.build_refmask()

        logger.info("Processing pileups...")
        for column in pileup_iter:
            scaffold = column.reference_name
            refpos = column.reference_pos

            for read in column.pileups:
                self._assess_read(call_data, scaffold, refpos, read)
                self._assess_alternative_locations(call_data, refpos, read)

        logger.info("Done.")
        logger.info("Analyzing coverage...")
        call_data.analyze_coverage()

        logger.info("Calling alleles...")
        call_data.call_alleles(self.min_pileup_qual, self.min_qual_frac)

        logger.info("Finding gaps...")
        call_data.find_gaps()
        logger.info("Done.")

        return call_data

    def _assess_read(self, call_data, scaffold, refpos, read):
        alignment = read.alignment

        # if this is a paired read, make sure the pairs are properly aligned
        if alignment.is_paired and not alignment.is_proper_pair:
            call_data.bad_read(scaffold, refpos)
            return

        # restrict ourselves to full-length alignments (not clipped)
        if alignment.query_alignment_length != alignment.query_length:
            # alignment is clipped
            call_data.bad_read(scaffold, refpos)
            return

        # check that inferred insert size is at least read length
        if alignment.is_paired:
            tlen = alignment.template_length
            if abs(tlen) < alignment.query_length:
                call_data.bad_read(scaffold, refpos)
                return

        # get base quality (note this is next base if deletion, but we won't
        # use that)
        pos = read.query_position_or_next
        qual = alignment.query_qualities[pos]
        if qual < self.min_qual:
            call_data.bad_read(scaffold, refpos)
            return

        # check for decent mapping quality
        mq = alignment.mapping_quality
        call_data.update_mapping_quality(scaffold, refpos, mq)

        # we keep track of otherwise good reads with low mapping quality;
        # that probably means this is a repeat
        if mq < self.min_mapping_quality:
            call_data.low_mapping_quality(scaffold, refpos)
            return

        if read.indel:
            if read.is_del:
                # deletion
                base = Allele.DEL
            else:
                # then it must be an insertion
                base = Allele.INS
        else:
            # base call must be real base (e.g., not N)
            base = Allele.from_str(alignment.query_sequence[pos])
            if not base:
                call_data.bad_read(scaffold, refpos)
                return

        # We're good! Update the pileup stats...
        call_data.good_read(scaffold, refpos, base, qual)

    def _assess_alternative_locations(self, call_data, refpos, read):
        """
        Assess alternative alignment locations of a read.

        We keep track of otherwise good reads with low mapping
        quality; that probably means this is a repeat.

        The code below updates "low mapping quality" counts for
        every location where this read maps equally well. This
        influences gap prediction.
        """

        alignment = read.alignment
        mq = alignment.mapping_quality
        if mq < self.min_mapping_quality:
            for scaffold, pos, rc in self._alternative_locations(alignment,
                                                                 refpos):
                call_data.low_mapping_quality(scaffold, pos)

    def _alternative_locations(self, read, loc):
        if read.has_tag("XA"):
            xa = read.get_tag("XA")
            nm = int(read.get_tag("NM"))
            offset = (read.reference_end - loc + 1 if read.is_reverse else
                      loc - read.reference_start)
            for aln in xa.split(';'):
                if not aln:
                    continue

                scaffold, pos, cigar, alt_nm = aln.split(',')

                if 'S' in cigar or 'H' in cigar:
                    # Clipped alignment, ignore
                    logger.debug("Ignoring clipped alternative alignment")
                    continue

                alt_nm = int(alt_nm)
                if alt_nm <= nm:
                    pos = int(pos)
                    rc = pos < 0

                    # Turn into a 0-based coordinate system
                    pos = abs(pos) - 1
                    coord = (pos + read.query_length - offset + 1 if rc
                             else pos + offset)

                    yield scaffold, coord, rc
