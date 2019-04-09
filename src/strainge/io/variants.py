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

import io
import csv
import logging
import itertools
from datetime import datetime

import h5py
import numpy
import vcf
from vcf.model import (_Record as VcfRecord, _Substitution as Substitution,
                       _SV as SV)

from strainge.utils import find_consecutive_groups
from strainge.variant_caller import Allele, VariantCallData

logger = logging.getLogger(__name__)

TSV_FIELDS = (
    ("name", "%s"), ("length", "%d"), ("coverage", "%.3f"),
    ("median", "%d"), ("callable", "%d"), ("callablePct", "%.3f"),
    ("confirmed", "%d"), ("confirmedPct", "%.3f"),
    ("snps", "%d"), ("snpPct", "%.3f"), ("pureSnps", "%d"),
    ("pureSnpPct", "%.3f"), ("lowmq", "%d"), ("lowmqPct", "%.3f"),
    ("high", "%d"), ("highPct", "%.3f"), ("gapCount", "%d"),
    ("gapLength", "%d")
)

VCF_TEMPLATE = """##fileformat=VCFv4.0
##fileDate={date}
##source=StrainGR
##reference={ref}
##INFO=<ID=DP,Number=1,Type=Integer,Description="Coverage depth">
##INFO=<ID=RQ,Number=1,Type=Integer,Description="Sum of reference base \
qualities">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Mean mapping quality">
##INFO=<ID=RF,Number=1,Type=Float,Description="Reference Fraction">
##INFO=<ID=AD,Number=R,Type=Integer,Description="Allele read depths">
##INFO=<ID=QS,Number=R,Type=Integer,Description="Sum of base qualities">
##INFO=<ID=BF,Number=R,Type=Float,Description="Quality weighted base \
frequencies">
##INFO=<ID=ST,Number=R,Type=Integer,Description="Which of the observed \
alleles have strong evidence">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DEL,Description="Deletion">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def generate_call_summary_tsv(call_data, output_file):
    """
    This function creates a TSV file containing summary variant calling
    statistics per scaffold as well for the complete genome as a whole.

    Parameters
    ----------
    call_data : VariantCallData
        `VariantCallData` object containing all variant calling statistics.
    output_file : file
        File-like object where the TSV will be written
    """

    writer = csv.writer(output_file, delimiter='\t', lineterminator='\n')

    # Write header
    writer.writerow(f[0] for f in TSV_FIELDS)

    # Write each row, `stats` is a dictionary with all relevant data per
    # scaffold
    for stats in call_data.summarize():
        writer.writerow(f[1] % stats.get(f[0], 0) for f in TSV_FIELDS)


def call_data_to_hdf5(call_data, output_file):
    """Write several summary statistics calculated by `VariantCallData`
    to a HDF5 file.

    Parameters
    ----------
    call_data : VariantCallData
    output_file : str
        Output file name
    """

    with h5py.File(output_file, 'w') as h5:
        for scaffold in call_data.scaffolds_data.values():
            scaffold_grp = h5.create_group(scaffold.name)

            scaffold_grp.create_dataset(
                "refmask", data=scaffold.refmask, compression=9)
            scaffold_grp.create_dataset(
                "alleles", data=scaffold.alleles, compression=9)
            scaffold_grp.create_dataset(
                "bad", data=scaffold.bad, compression=9)
            scaffold_grp.create_dataset(
                "lowmq_count", data=scaffold.lowmq_count, compression=9)
            scaffold_grp.create_dataset(
                "mq_sum", data=scaffold.mq_sum, compression=9)
            scaffold_grp.create_dataset(
                "strong", data=scaffold.strong, compression=9)
            scaffold_grp.create_dataset(
                "weak", data=scaffold.weak, compression=9)
            scaffold_grp.create_dataset(
                "coverage", data=scaffold.coverage, compression=9)
            scaffold_grp.create_dataset(
                "high_coverage", data=scaffold.high_coverage, compression=9)

            scaffold_grp.attrs["mean_coverage"] = scaffold.mean_coverage
            scaffold_grp.attrs["median_coverage"] = scaffold.median_coverage
            scaffold_grp.attrs["coverage_cutoff"] = scaffold.coverage_cutoff

        h5.attrs['type'] = "VariantCallData"
        h5.attrs['min_gap_size'] = call_data.min_gap_size
        h5.attrs['mean_coverage'] = call_data.mean_coverage
        h5.attrs['median_coverage'] = call_data.median_coverage
        h5.attrs['reference_fasta'] = call_data.reference_fasta


def call_data_from_hdf5(hdf5_file):
    """
    Create a `CallStatsCollector` by loading the relevant data from an
    earlier created HDF5 file.

    Parameters
    ----------
    hdf5_file : str
        HDF5 filename

    Returns
    -------
    CallStatsCollector
    """
    datasets = {"refmask", "alleles", "bad", "lowmq_count",
                "mq_sum", "strong", "weak", "coverage", "high_coverage"}

    # These datasets have pre-allocated numpy arrays in `VariantCallData`
    read_direct = {"alleles", "bad", "lowmq_count", "mq_sum"}

    with h5py.File(hdf5_file, 'r') as hdf5:
        if 'type' not in hdf5.attrs:
            raise IOError(f"The HDF5 file {hdf5_file} does not contain"
                          f"`VariantCallData`.")

        min_gap = hdf5.attrs['min_gap_size']

        scaffolds = {}
        for scaffold in hdf5:
            scaffolds[scaffold] = len(hdf5[scaffold]['refmask'])

        call_data = VariantCallData(scaffolds, min_gap)

        for scaffold_name in hdf5:
            scaffold = call_data.scaffolds_data[scaffold_name]

            for dataset_name in hdf5[scaffold_name]:
                if dataset_name not in datasets:
                    continue

                logger.info("Loading dataset '%s' for scaffold '%s' from "
                            "HDF5...", dataset_name, scaffold_name)

                if dataset_name in read_direct:
                    target = getattr(scaffold, dataset_name)
                    hdf5[scaffold_name][dataset_name].read_direct(target)
                else:
                    arr = numpy.array(hdf5[scaffold_name][dataset_name])
                    setattr(scaffold, dataset_name, arr)

            scaffold.mean_coverage = hdf5[scaffold_name].attrs['mean_coverage']
            scaffold.median_coverage = hdf5[scaffold_name].attrs[
                'median_coverage']
            scaffold.coverage_cutoff = hdf5[scaffold_name].attrs[
                'coverage_cutoff']

        call_data.mean_coverage = hdf5.attrs['mean_coverage']
        call_data.median_coverage = hdf5.attrs['median_coverage']
        call_data.min_gap_size = hdf5.attrs['min_gap_size']

        if 'reference_fasta' in hdf5.attrs:
            call_data.reference_fasta = hdf5.attrs['reference_fasta']
        else:
            call_data.reference_fasta = ""

        logger.info("Mean coverage (across all scaffolds): %.2f",
                    call_data.mean_coverage)
        logger.info("Median coverage (across all scaffolds): %.2f",
                    call_data.median_coverage)
        logger.info("Minimum gap size: %.2f", call_data.min_gap_size)

        # Reconstruct gaps again
        call_data.find_gaps()

    return call_data


def boolean_array_to_bedfile(array, output_file, scaffold_name,
                             min_feature_size=1):
    """Convert a boolean numpy array to a BED file, which can be visualised in
    genome browsers. This function searches for groups of consecutive 1's
    (True's), and each such group is written to the file as a single feature.

    Parameters
    ----------
    array : ndarray
        The one dimensional boolean array to process
    output_file : file-like object
        Output file object where the data will be written to.
    scaffold_name : str
        The name of the scaffold the array corresponds to.
    min_feature_size : int
        Only output the feature if it larger than the given size. Defaults
        to 1.
    """

    writer = csv.writer(output_file, delimiter='\t', lineterminator='\n')

    for group, start, end, length in find_consecutive_groups(
            array, min_feature_size):
        if not numpy.all(group):
            continue

        writer.writerow((scaffold_name, start, end))


def array_to_bedgraph(array, output_file, scaffold_name):
    """
    Export a numpy array to BedGraph format, which can be visualized in
    genome viewers.

    Parameters
    ----------
    array : ndarray
    output_file : file-like object
    scaffold_name : str
    """

    writer = csv.writer(output_file, delimiter='\t', lineterminator='\n')

    for group, start, end, length in find_consecutive_groups(array):
        writer.writerow([scaffold_name, start, end, group[0]])


def vcf_records_for_scaffold(scaffold, verboseness=0):
    """
    Yield VCF records for all positions in the scaffold, or only the
    positions with something else than the reference.

    Parameters
    ----------
    scaffold : strainge.variant_caller.ScaffoldCallData
        Variant call data for this scaffold
    verboseness : int
        Determines the verboseness of the VCF. Higher value will result in
        more entries in the VCF. Accepted values:

        0: Only output StrainGR strong SNPs
        1: Output strong and weak SNPs
        2: Output an entry for every position in the genome, even if nothing
           else but the reference is observed.

    """

    logger.info("Generate VCF records for scaffold %s", scaffold.name)
    positions = numpy.arange(scaffold.length)

    if verboseness < 2:
        # Find positions with something else than the reference base
        # Remove bit corresponding to reference base
        if verboseness == 0:
            alt_mask = scaffold.strong & ~scaffold.refmask
        else:
            alt_mask = scaffold.weak & ~scaffold.refmask

        # Only keep those with something else than the reference
        ix = alt_mask > 0
        positions = positions[ix]

    for pos in positions:
        alts = []
        alts_bit = []
        strong = set()

        # Check for alternative alleles (i.e. not the reference base)
        for allele in Allele:
            if (allele & scaffold.weak[pos] and
                    not allele & scaffold.refmask[pos]):
                alts_bit.append(allele)

                if allele & (Allele.INS | Allele.DEL):
                    alts.append(SV(str(allele)))
                else:
                    alts.append(Substitution(str(allele)))

            # We don't check if this is the refbase or not because we also want
            # to include information whether the reference is strongly
            # confirmed or not.
            if allele & scaffold.strong[pos]:
                strong.add(allele)

        ref_plus_alts = [scaffold.refmask[pos]] + alts_bit
        allele_counts = [scaffold.allele_count(pos, allele)
                         for allele in ref_plus_alts]
        allele_quals = [scaffold.allele_qual(pos, allele)
                        for allele in ref_plus_alts]

        yield VcfRecord(
            CHROM=scaffold.name,
            POS=pos+1,  # 1-based coordinate system
            ID=".",
            REF=str(Allele(scaffold.refmask[pos]))[:1],
            ALT=alts,
            QUAL=".",
            FILTER="PASS",
            INFO={
                'DP': scaffold.depth(pos),
                'MQ': int(round(scaffold.mean_mq(pos))),
                'RQ': scaffold.ref_qual(pos),
                'RF': round(scaffold.ref_fraction(pos), 3),
                'AD': allele_counts,
                'QS': allele_quals,
                'ST': [int(b in strong) for b in ref_plus_alts]
            },
            FORMAT="",
            sample_indexes=None
        )


def write_vcf(call_data, output_file, verboseness=0):
    """
    Write out variant calling data to a VCF file. By default only writes
    entries where something else than the reference has been observed.

    Parameters
    ----------
    call_data : VariantCallData
    output_file : file
    verboseness : int
        Determines the verboseness of the VCF. Higher value will result in
        more entries in the VCF. Accepted values:

        0: Only output StrainGR strong SNPs
        1: Output strong and weak SNPs
        2: Output an entry for every position in the genome, even if nothing
           else but the reference is observed.
    """

    vcf_template = vcf.Reader(io.StringIO(VCF_TEMPLATE.format(
        date=datetime.now(),
        ref=call_data.reference_fasta
    )))

    vcf_writer = vcf.Writer(output_file, vcf_template)

    record_iter = itertools.chain.from_iterable(
        vcf_records_for_scaffold(scaffold, verboseness)
        for scaffold in call_data.scaffolds_data.values()
    )

    for record in record_iter:
        vcf_writer.write_record(record)
