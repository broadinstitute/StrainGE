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
import json
from typing import List, Iterable  # noqa
from pathlib import Path
from collections import defaultdict

from Bio import SeqIO

from strainge.kmertools import open_seq_file


def parse_straingst(result_file, return_sample_stats=False):
    """Parse StrainGST output file and return the strains present in a sample
    along with all metrics.

    Returns
    -------
    Iterable[dict]
    """

    # Ignore comments
    result_file = (line for line in result_file if not line.startswith('#'))

    # Collect sample statistics (first two lines)
    sample_stats = [
        next(result_file),
        next(result_file)
    ]

    if return_sample_stats:
        sample_stats = next(csv.DictReader(sample_stats, delimiter='\t'))

        # Return sample statistics
        yield sample_stats

    # Return each strain found with its statistics
    yield from csv.DictReader(result_file, delimiter='\t')


def ref_concat_with_metadata(refs, concat_out, metadata_out=None):
    """Concatenate the given list of references, and keep track which contigs
    belong to which reference. This will be written to `metadata_out`.

    Parameters
    ----------
    refs : List[str]
    concat_out : file
    metadata_out : file
    """

    contigs = defaultdict(list)
    to_write = []
    for ref in refs:
        ref_path = Path(ref)

        for scaffold in open_seq_file(ref):
            contigs[ref_path.stem].append(scaffold.name)
            to_write.append(scaffold)

    SeqIO.write(to_write, concat_out, "fasta")

    if metadata_out:
        json.dump(contigs, metadata_out)
