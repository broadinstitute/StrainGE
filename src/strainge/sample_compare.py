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

import numpy
from intervaltree import IntervalTree

from strainge.variant_caller import count_ts_tv, scale_min_gap_size
from strainge.utils import pct

logger = logging.getLogger(__name__)


class SampleComparison:
    """
    This class compares variant calls in two different samples, and gives
    statistics on how similar the strains are.
    """

    def __init__(self, call_data1, call_data2):
        """
        Compare the variant call data from two samples.

        Parameters
        ----------
        call_data1 : strainge.variant_caller.VariantCallData
            Variant call data of sample 1
        call_data2 : strainge.variant_caller.VariantCallData
            Variant call data of sample 2
        """

        scaffolds_common = (call_data1.scaffolds_data.keys() &
                            call_data2.scaffolds_data.keys())

        logger.info("Scaffolds in common: %s", scaffolds_common)

        self.metrics = {}

        self.sample1 = call_data1
        self.sample2 = call_data2

        for scaffold in scaffolds_common:
            scaffold_a = call_data1.scaffolds_data[scaffold]
            scaffold_b = call_data2.scaffolds_data[scaffold]

            assert scaffold_a.length == scaffold_b.length

            self.metrics[scaffold] = self._do_compare(scaffold_a, scaffold_b)
            self.metrics[scaffold].update(self.compare_gaps(scaffold_a,
                                                            scaffold_b))

    def _do_compare(self, a, b):
        """

        Parameters
        ----------
        a : strainge.variant_caller.ScaffoldCallData
        b : strainge.variant_caller.ScaffoldCallData

        Returns
        -------
        dict
        """
        # common locations where both have a call
        common, common_cnt, common_pct = self.compare_thing(
            numpy.ones_like(a.refmask), numpy.logical_and(a.strong, b.strong))

        # locations where both have only a single allele called
        single_a = (a.strong & (a.strong - 1)) == 0
        single_b = (b.strong & (b.strong - 1)) == 0
        singles, single_cnt, single_pct = self.compare_thing(
            common, single_a & single_b)

        # locations where both have the same single allele called
        single_agree, single_agree_cnt, single_agree_pct = self.compare_thing(
            singles, a.strong == b.strong)

        _, multi_cnt, multi_pct = self.compare_thing(
            common, ~single_a | ~single_b)

        # Common locations where they share at least one allele (other alleles
        # may be present)
        _, shared_alleles_cnt, shared_alleles_pct = self.compare_thing(
            common, (a.strong & b.strong) > 0)

        # common locations where either has a variant from reference
        variants, variant_cnt, variant_pct = self.compare_thing(
            common, ((a.strong | b.strong) & ~a.refmask) > 0)

        # variant locations where both have a shared variant
        common_var, common_var_cnt, common_var_pct = self.compare_thing(
            variants, (a.strong & b.strong) > 0)

        # variant locations where both agree
        var_agree, var_agree_cnt, var_agree_pct = self.compare_thing(
            variants, a.strong == b.strong)

        # variant in a but not b
        a_not_b, a_not_b_cnt, a_not_b_pct = self.compare_thing(
            variants, (a.strong & ~b.strong & ~a.refmask) > 0)

        # variant in a but not b weakly
        a_not_bweak, a_not_bweak_cnt, a_not_bweak_pct = self.compare_thing(
            variants, (a.strong & ~b.weak & ~a.refmask) > 0)

        # variant in b not a
        b_not_a, b_not_a_cnt, b_not_a_pct = self.compare_thing(
            variants, (b.strong & ~a.strong & ~a.refmask) > 0)

        # variant in b not a weakly
        b_not_aweak, b_not_aweak_cnt, b_not_aweak_pct = self.compare_thing(
            variants, (b.strong & ~a.weak & ~a.refmask) > 0)

        # Count transitions/transversions
        disagree, disagree_cnt, disagree_pct = self.compare_thing(
            singles, a.strong != b.strong)

        transitions, transversions = count_ts_tv(a.strong[disagree],
                                                 b.strong[disagree])

        transitions_pct = transitions / single_cnt if single_cnt else 0.0
        transversions_pct = transitions / single_cnt if single_cnt else 0.0

        return {
            "ref": a.ref_name if a.ref_name else "na",
            "length": a.length,
            "common": common_cnt,
            "commonPct": common_pct,
            "single": single_cnt,
            "singlePct": single_pct,
            "singleAgree": single_agree_cnt,
            "singleAgreePct": single_agree_pct,
            "multi": multi_cnt,
            "multiPct": multi_pct,
            "sharedAlleles": shared_alleles_cnt,
            "sharedAllelesPct": shared_alleles_pct,
            "variants": variant_cnt,
            "variantPct": variant_pct,
            "commonVariant": common_var_cnt,
            "commonVariantPct": common_var_pct,
            "variantExact": var_agree_cnt,
            "variantExactPct": var_agree_pct,
            "AnotB": a_not_b_cnt,
            "AnotBpct": a_not_b_pct,
            "AnotBweak": a_not_bweak_cnt,
            "AnotBweakPct": a_not_bweak_pct,
            "BnotA": b_not_a_cnt,
            "BnotApct": b_not_a_pct,
            "BnotAweak": b_not_aweak_cnt,
            "BnotAweakPct": b_not_aweak_pct,
            "transitions": transitions,
            "tsPct": transitions_pct * 100,
            "transversions": transversions,
            "tvPct": transversions_pct * 100
        }

    def compare_thing(self, common, thing):
        """
        Computes occurrence of a condition within a set, and returns those
        stats.

        :param common: flag for locations to consider
        :param thing: condition we're looking for in common

        :return: array where command and condition are true, count of that,
                 and percentage with respect to common
        """
        common_cnt = numpy.count_nonzero(common)
        common_things = numpy.logical_and(common, thing)
        common_things_cnt = numpy.count_nonzero(common_things)
        percent = pct(common_things_cnt, common_cnt)

        return common_things, common_things_cnt, percent

    def compare_gaps(self, a, b):
        """
        More stats, this time about uncovered regions in common.

        Parameters
        ----------
        a : strainge.variant_caller.ScaffoldCallData
        b : strainge.variant_caller.ScaffoldCallData
        """

        a_length = sum(g.length for g in a.gaps)
        b_length = sum(g.length for g in b.gaps)

        gap_tree_a = IntervalTree.from_tuples(
            (g.start, g.end, g) for g in a.gaps
        )
        gap_tree_b = IntervalTree.from_tuples(
            (g.start, g.end, g) for g in b.gaps
        )

        a_shared = [g for g in a.gaps if gap_tree_b[g.start:g.end]]
        b_shared = [g for g in b.gaps if gap_tree_a[g.start:g.end]]

        a_shared_length = sum(g.length for g in a_shared)
        b_shared_length = sum(g.length for g in b_shared)

        # Calculate gap similarity between two samples using jaccard
        gaps_a = numpy.zeros((a.length,), dtype=bool)
        gaps_b = numpy.zeros((b.length,), dtype=bool)

        # For jaccard similarity, only take into account gaps that can be
        # detected in both samples. That means, that if one sample has a
        # higher min_gap_size because it's a low abundance strain,
        # ignore gaps in the other sample smaller than that size.
        scaled_min_gap_a = scale_min_gap_size(self.sample1.min_gap_size,
                                              a.mean_coverage)
        scaled_min_gap_b = scale_min_gap_size(self.sample2.min_gap_size,
                                              b.mean_coverage)
        min_gap_size = max(scaled_min_gap_a, scaled_min_gap_b)

        for g in a.gaps:
            if (g.end - g.start) >= min_gap_size:
                gaps_a[g.start:g.end] = True

        for g in b.gaps:
            if (g.end - g.start) >= min_gap_size:
                gaps_b[g.start:g.end] = True

        total_non_gap = (~gaps_a | ~gaps_b).sum()
        if total_non_gap > 0:
            jaccard = (~gaps_a & ~gaps_b).sum() / total_non_gap
        else:
            jaccard = 1.0

        return {
            "Agaps": a_length,
            "AsharedGaps": a_shared_length,
            "AgapPct": pct(a_shared_length, a_length),
            "Bgaps": b_length,
            "BsharedGaps": b_shared_length,
            "BgapPct": pct(b_shared_length, b_length),
            "gapJaccardSim": jaccard,
        }
