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

from strainge import kmerizer


def jaccard(kmers1, kmers2):
    """Computes jaccard similarity. Returns numerator and denominator
    separately."""
    intersection = kmerizer.count_common(kmers1, kmers2)
    return intersection / (kmers1.size + kmers2.size - intersection)


def minsize(kmers1, kmers2):
    intersection = kmerizer.count_common(kmers1, kmers2)
    return intersection / min(kmers1.size, kmers2.size)


def meansize(kmers1, kmers2):
    intersection = kmerizer.count_common(kmers1, kmers2)
    return intersection / ((kmers1.size + kmers2.size) / 2)


def maxsize(kmers1, kmers2):
    intersection = kmerizer.count_common(kmers1, kmers2)
    return intersection / max(kmers1.size, kmers2.size)


def subset(kmers1, kmers2):
    """Calculate the fraction of k-mers in k-merset 1 that are also in k-merset
    2, useful to check whether k-merset 1 is a subset of another."""
    intersection = kmerizer.count_common(kmers1, kmers2)
    return intersection / kmers1.size


def reference(kmers1, kmers2):
    """Assume k-merset 2 is the k-merset of a reference genome."""
    intersection = kmerizer.count_common(kmers1, kmers2)
    return intersection / kmers2.size


def similarity_score(kmers1, kmers2, scoring="jaccard"):
    """Compute how similar two k-mer sets are using various scoring methods.
    Returns numerator and denominator separately."""

    if scoring not in SCORING_METHODS:
        raise ValueError("Invalid scoring method '{}'".format(scoring))

    return SCORING_METHODS[scoring](kmers1, kmers2)


def ani(k, j):
    """Estimate average nucleotide identity from Jaccard distance between
    two k-mer sets. Also known as mash [1] distance.

    [1]: Mash: fast genome and metagenome distance estimation using MinHash.
         Ondov BD, Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S,
         Phillippy AM. Genome Biol. 2016 Jun 20;17(1):132.
         doi: 10.1186/s13059-016-0997-x.
    """

    distance = -math.log(2 * j / (1.0 + j)) / k if j > 0 else 1

    # Transform to ANI
    return 1 - distance


SCORING_FUNCS = (jaccard, minsize, meansize, maxsize, subset, reference)
SCORING_METHODS = {func.__name__: func for func in SCORING_FUNCS}
