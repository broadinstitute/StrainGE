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
import bz2
import gzip
from typing import List, Iterable  # noqa
from pathlib import Path
from contextlib import contextmanager
from itertools import chain


@contextmanager
def open_compressed(filename):
    if not isinstance(filename, Path):
        filename = Path(filename)

    if filename.suffix == ".gz":
        f = gzip.open(filename, "rt")
    elif filename.suffix == ".bz2":
        f = bz2.open(filename, "rt")
    else:
        f = open(filename)

    yield f

    f.close()


def read_fastq(fp):
    """Heng Li's fast FASTQ reader."""

    last = None
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break

        if not last:
            break

        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])

        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break


def parse_straingst(result_file, return_sample_stats=False):
    """Parse StrainGST output file and return the strains present in a sample
    along with all metrics.

    Returns
    -------
    Iterable[dict]
    """

    # Ignore comments
    result_file = (line for line in result_file if not line.startswith('#'))
    first_line = next(result_file)

    old_style_straingst = False
    if first_line.startswith("sample"):
        old_style_straingst = True

    # Collect sample statistics (first two lines)
    if old_style_straingst:
        sample_stats = [
            first_line,
            next(result_file)
        ]
    else:
        sample_stats = []

    if sample_stats and return_sample_stats:
        sample_stats = next(csv.DictReader(sample_stats, delimiter='\t'))

        # Return sample statistics
        yield sample_stats
    elif not sample_stats and return_sample_stats:
        raise ValueError("Trying to read sample statistics from a new style (v1.3) StrainGST output!")

    # Return each strain found with its statistics
    if old_style_straingst:
        yield from csv.DictReader(result_file, delimiter='\t')
    else:
        yield from csv.DictReader(chain([first_line], result_file), delimiter='\t')
