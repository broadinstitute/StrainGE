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

import csv

from strainge.variant_caller import Allele

COMPARE_TSV_FIELDS = (
    ("scaffold", "%s"), ("common", "%d"),
    ("commonPct", "%.3f"), ("single", "%d"),
    ("singlePct", "%.3f"), ("singleAgree", "%d"), ("singleAgreePct", "%.3f"),
    ("variants", "%d"), ("variantPct", "%.3f"), ("commonVariant", "%d"),
    ("commonVariantPct", "%.3f"), ("variantAgree", "%d"),
    ("variantAgreePct", "%.3f"),
    ("transitions", "%d"), ("tsPct", "%.3f"),
    ("transversions", "%d"), ("tvPct", "%.3f"),
    ("AnotB", "%d"), ("AnotBpct", "%.3f"),
    ("AnotBweak", "%d"), ("AnotBweakPct", "%.3f"), ("BnotA", "%d"),
    ("BnotApct", "%.3f"), ("BnotAweak", "%d"), ("BnotAweakPct", "%.3f"),
    ("Agaps", "%d"), ("AsharedGaps", "%d"), ("AgapPct", "%.3f"),
    ("Bgaps", "%d"), ("BsharedGaps", "%d"), ("BgapPct", "%.3f"),
    ("gapJaccardSim", "%.3f"),
)


def generate_compare_summary_tsv(sample_comparison, output_file):
    """
    Generate a TSV summarizing the results of a sample comparison

    Parameters
    ----------
    sample_comparison : strainge.sample_compare.SampleComparison
    output_file : file object
    """

    writer = csv.writer(output_file, delimiter='\t', lineterminator='\n')

    writer.writerow(f[0] for f in COMPARE_TSV_FIELDS)

    for scaffold, metrics in sample_comparison.metrics.items():
        metrics['scaffold'] = scaffold
        writer.writerow([f[1] % metrics.get(f[0], 0)
                         for f in COMPARE_TSV_FIELDS])


def generate_compare_details_tsv(f, a, b, verbose=False):
    """
    Generate a TSV describing variant call differences between the two given
    samples, down to the single nucleotide level.

    Parameters
    ----------
    f : file
        Output file object
    a : strainge.variant_caller.VariantCallData
    b : strainge.variant_caller.VariantCallData
    """

    writer = csv.writer(f, delimiter='\t', lineterminator='\n')

    for scaffoldA, scaffoldB in zip(a.scaffolds_data.values(),
                                    b.scaffolds_data.values()):
        assert scaffoldA.name == scaffoldB.name
        assert scaffoldA.length == scaffoldB.length

        for pos in range(scaffoldA.length):
            strongA = scaffoldA.strong[pos]
            strongB = scaffoldB.strong[pos]

            if strongA and strongB and (verbose or strongA != strongB):
                ref = Allele(scaffoldA.refmask[pos])
                alleles_a = _strong_weak_string(strongA, scaffoldA.weak[pos])
                alleles_b = _strong_weak_string(strongB, scaffoldB.weak[pos])
                writer.writerow([scaffoldA.name, pos, str(ref), alleles_a,
                                 alleles_b])


def _strong_weak_string(strong, weak):
    """
    Returns a string describing strong and weak allele calls.

    Strong are weak calls are separated by a /, and both strong and weak can
    have multiple alleles. Example:

        A/C,T

    A is a strong call, while C and T are weak calls.
    """

    strong = Allele(int(strong))
    weak = Allele(int(weak & ~strong))

    if weak:
        return f"{str(strong)}/{str(weak)}"
    else:
        return str(strong)
