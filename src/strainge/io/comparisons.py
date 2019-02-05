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

COMPARE_TSV_FIELDS = (
    ("scaffold", "%s"), ("common", "%d"),
    ("commonPct", "%.2f"), ("single", "%d"),
    ("singlePct", "%.2f"), ("singleAgree", "%d"), ("singleAgreePct", "%.2f"),
    ("variants", "%d"), ("variantPct", "%.2f"), ("commonVariant", "%d"),
    ("commonVariantPct", "%.2f"), ("variantAgree", "%d"),
    ("variantAgreePct", "%.2f"), ("AnotB", "%d"), ("AnotBpct", "%.2f"),
    ("AnotBweak", "%d"), ("AnotBweakPct", "%.2f"), ("BnotA", "%d"),
    ("BnotApct", "%.2f"), ("BnotAweak", "%d"), ("BnotAweakPct", "%.2f"),
    ("Agaps", "%d"), ("AsharedGaps", "%d"), ("AgapPct", "%.2f"),
    ("Bgaps", "%d"), ("BsharedGaps", "%d"), ("BgapPct", "%.2f")
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
        writer.writerow([scaffold] + [f[1] % metrics.get(f[0], 0)
                                      for f in COMPARE_TSV_FIELDS])


