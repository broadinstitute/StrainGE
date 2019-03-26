"""
StrainGE utility functions
==========================
"""

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
from collections import namedtuple

import numpy


Group = namedtuple('Group', 'data start end length')


def parse_num_suffix(num):
    """
    Parse a string containing a number with a possible suffix like 'M' or 'G',
    and multiply accordingly.

    >>> parse_num_suffix('10M')
    10000000
    >>> parse_num_suffix('5G')
    5000000000
    >>> parse_num_suffix('3k')
    3000
    >>> parse_num_suffix('500')
    500

    Parameters
    ----------
    num : str
        The number with possible suffix to parse

    Returns
    -------
    int
        An integer multiplied accordingly to the suffix
    """

    if not num:
        return None

    suffixes = {
        'G': 1000000000,
        'M': 1000000,
        'K': 1000
    }

    if not num[-1].isalpha():
        return int(num)

    suffix = num[-1].upper()
    if suffix not in suffixes:
        raise ValueError(
            "'{}' is not a valid number. Supported suffixes: {}".format(
                num, ", ".join(iter(suffixes.keys()))
            ))

    return int(num[:-1]) * suffixes[suffix]


def pct(numerator, denominator, precision=None):
    """Makes into a percent, avoiding division by zero"""

    if numerator > 0 and denominator > 0:
        value = (100.0 * numerator) / denominator
    else:
        value = 0.0

    if precision is not None:
        value = round(value, precision)

    return value


def lander_waterman(coverage):
    """
    Use Lander-Waterman relation to compute expect fraction of genome
    covered given mean coverage.

    Parameters
    ----------
    coverage : float
        Mean coverage

    Returns
    -------
    Expected fraction of genome covered
    """
    return 1.0 - math.exp(-coverage)


def find_consecutive_groups(array, min_size=1):
    """
    Find consecutive groups numbers in a numpy boolean array.
    Adjusted from https://stackoverflow.com/questions/7352684/

    Parameters
    ----------
    array : array-like
        Numpy 1-dimensional boolean array to process
    min_size
        Minimum size of the group to be reported

    Yields
    ------
    Group
        `Group` is a named tuple with attributes `data`, `start`, `end` and
        `size`.
    """

    # Positions where the number changes, i.e. derivative != 0
    change = numpy.diff(array) != 0

    # split at changing positions
    groups = numpy.split(array, numpy.where(change)[0] + 1)

    cur_pos = 0
    for group in groups:
        start = cur_pos
        length = len(group)
        cur_pos += length

        if length >= min_size:
            yield Group(group, start, cur_pos, length)

