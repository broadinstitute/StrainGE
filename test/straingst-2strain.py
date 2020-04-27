#!/usr/bin/env python
#  Copyright (c) 2016-2020, Broad Institute, Inc. All rights reserved.
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
import random


def cov_combinations():
    coverages = ('10x', '1x', '0.5x', '0.1x')
    return [(coverages[i], coverages[j]) for i in range(len(coverages)) for j in range(i, len(coverages))]


def gen_pairs():
    indb = list(range(1, 101))
    outdb = list(range(101, 201))

    random.shuffle(indb)
    random.shuffle(outdb)

    pairs = []
    # first, we take the first 25 pairs from same list
    for i in range(0, 50, 2):
        pairs.append((indb[i], indb[i+1]))
        pairs.append((outdb[i], outdb[i+1]))

    # then pick the next 50 for cross-db, alternating whether the first comes from indb or outdb
    for i in range(50, 100):
        if i % 2 == 0:
            pairs.append((indb[i], outdb[i]))
        else:
            pairs.append((outdb[i], indb[i]))

    return pairs


def sample_name(sample, cov):
    return f"sample{sample}-{cov}"


def gen_all():
    for cc in cov_combinations():
        outdir = f"{cc[0]}-{cc[1]}"
        print(f"mkdir {outdir}")
        for pair in gen_pairs():
            sample1 = sample_name(pair[0], cc[0])
            sample2 = sample_name(pair[1], cc[1])
            print(f"qs -p 2 straingst kmermerge -f -o {outdir}/{sample1}-{sample2}.hdf5 {cc[0]}/{sample1}-bg.hdf5 {cc[1]}/{sample2}.hdf5")

gen_all()

