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
import logging
import itertools
from datetime import datetime

import h5py
import numpy
import pysam

from strainge.utils import find_consecutive_groups
from strainge.variant_caller import Allele, VariantCallData

logger = logging.getLogger(__name__)

TSV_FIELDS = (
    ("ref", "%s"), ("name", "%s"), ("length", "%d"), ("coverage", "%.3f"),
    ("uReads", "%d"), ("abundance", "%.3f"),
    ("median", "%d"), ("callable", "%d"), ("callablePct", "%.3f"),
    ("confirmed", "%d"), ("confirmedPct", "%.3f"),
    ("snps", "%d"), ("snpPct", "%.3f"), ("multi", "%d"),
    ("multiPct", "%.3f"), ("lowmq", "%d"), ("lowmqPct", "%.3f"),
    ("high", "%d"), ("highPct", "%.3f"), ("gapCount", "%d"),
    ("gapLength", "%d")
)

VCF_TEMPLATE = """##fileformat=VCFv4.0
##fileDate={date}
##source=StrainGR
##reference={ref}
{contig_lengths}
##INFO=<ID=DP,Number=1,Type=Integer,Description="Coverage depth">
##INFO=<ID=RQ,Number=1,Type=Integer,Description="Sum of reference base qualities">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Mean mapping quality">
##INFO=<ID=RF,Number=1,Type=Float,Description="Reference Fraction">
##INFO=<ID=AD,Number=R,Type=Integer,Description="Allele read depths">
##INFO=<ID=QS,Number=R,Type=Integer,Description="Sum of base qualities">
##INFO=<ID=BF,Number=R,Type=Float,Description="Quality weighted base frequencies">
##INFO=<ID=ST,Number=R,Type=Integer,Description="Which of the observed alleles have strong evidence">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DEL,Description="Deletion">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"""


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
            scaffold_grp.attrs["read_count"] = scaffold.read_count
            scaffold_grp.attrs["repetitiveness"] = scaffold.repetitiveness
            scaffold_grp.attrs["ref_name"] = scaffold.ref_name

        h5.attrs['type'] = "VariantCallData"
        h5.attrs['min_gap_size'] = call_data.min_gap_size
        h5.attrs['mean_coverage'] = call_data.mean_coverage
        h5.attrs['median_coverage'] = call_data.median_coverage
        h5.attrs['reference_fasta'] = call_data.reference_fasta

        h5.attrs['total_reads'] = call_data.total_reads
        h5.attrs['passing_reads'] = call_data.passing_reads
        h5.attrs['lowmq_reads'] = call_data.lowmq_reads


def call_data_from_hdf5(hdf5_file, new_min_gap=None) -> VariantCallData:
    """
    Create a `CallStatsCollector` by loading the relevant data from an
    earlier created HDF5 file.

    Parameters
    ----------
    hdf5_file : str
        HDF5 filename
    new_min_gap : int
        Optionally set a new minimum gap size

    Returns
    -------
    VariantCallData
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
        if new_min_gap is not None:
            logger.info("Using new minimum gap size %d (original: %d)",
                        new_min_gap, min_gap)

            min_gap = new_min_gap

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

            if 'read_count' not in hdf5[scaffold_name].attrs:
                logger.warning("This is an old StrainGR hdf5 file. Strain "
                               "abundance information not available.")

            scaffold.read_count = hdf5[scaffold_name].attrs.get(
                'read_count', 0)

            if 'repetitiveness' not in hdf5[scaffold_name].attrs:
                logger.warning("This is an old StrainGR hdf5 file. Scaffold "
                               "repetitiveness information not available.")

            scaffold.repetitiveness = hdf5[scaffold_name].attrs.get(
                'repetitiveness', 0.0)

            scaffold.ref_name = hdf5[scaffold_name].attrs.get(
                'ref_name', "na")

        call_data.mean_coverage = hdf5.attrs['mean_coverage']
        call_data.median_coverage = hdf5.attrs['median_coverage']

        if 'total_reads' not in hdf5.attrs:
            logger.warning("This is an old StrainGR hdf5 file. Strain "
                           "abundance information not available.")

        call_data.total_reads = hdf5.attrs.get('total_reads', 0)

        if 'passing_reads' not in hdf5.attrs:
            logger.warning("This is an old StrainGR hdf5 file. Strain "
                           "abundance information not available.")

        call_data.passing_reads = hdf5.attrs.get('passing_reads', 0)

        if 'lowmq_reads' not in hdf5.attrs:
            logger.warning("This is an old StrainGR hdf5 file. Strain "
                           "abundance information not available.")

        call_data.lowmq_reads = hdf5.attrs.get('lowmq_reads', 0)

        if 'reference_fasta' in hdf5.attrs:
            call_data.reference_fasta = hdf5.attrs['reference_fasta']
        else:
            call_data.reference_fasta = ""

        logger.info("Mean coverage (across all scaffolds): %.2f",
                    call_data.mean_coverage)
        logger.info("Median coverage (across all scaffolds): %.2f",
                    call_data.median_coverage)

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


def array_to_wig(array, output_file, scaffold_name):
    """
    Write the values of a numpy array to WIG format, which is viewable in
    IGV. WIG is better suitable for continuous data.

    Parameters
    ----------
    array : array (1D)
    output_file : str or file-like object
        If the given filename contains ".gz" it will automatically compress
        the file.
    scaffold_name : str
    """

    numpy.savetxt(output_file, array.T, fmt="%g",
                  header=f"fixedStep chrom={scaffold_name}", comments="",
                  encoding="utf-8")


def vcf_records_for_scaffold(writer, scaffold, verboseness=0):
    """
    Yield VCF records for all positions in the scaffold, or only the
    positions with something else than the reference.

    Parameters
    ----------
    writer : pysam.VariantFile
        PySam output VCF object
    scaffold : strainge.variant_caller.ScaffoldCallData
        Variant call data for this scaffold
    verboseness : int
        Determines the verboseness of the VCF. Higher value will result in
        more entries in the VCF. Accepted values:

        0: Only output StrainGR strong SNPs
        1: Output strong and weak SNPs

    """

    logger.info("Generate VCF records for scaffold %s", scaffold.name)
    positions = numpy.arange(scaffold.length)

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
            called_alleles = (scaffold.strong[pos] if verboseness == 0 else
                              scaffold.weak[pos])
            if (allele & called_alleles and
                    not allele & scaffold.refmask[pos]):
                alts_bit.append(allele)

                if allele & (Allele.INS | Allele.DEL):
                    alts.append(str(allele))
                else:
                    alts.append(str(allele))

            # We don't check if this is the refbase or not because we also want
            # to include information whether the reference is strongly
            # confirmed or not.
            if allele & scaffold.strong[pos]:
                strong.add(allele)

        ref = scaffold.refmask[pos]
        ref_plus_alts = [Allele(ref)] + alts_bit

        allele_counts = [
            int(scaffold.allele_count(pos, allele)) if allele else 0
            for allele in ref_plus_alts
        ]
        allele_quals = [
            int(scaffold.allele_qual(pos, allele)) if allele else 0
            for allele in ref_plus_alts
        ]

        sum_quals = sum(allele_quals)
        if sum_quals:
            weighted_base_freqs = [v / sum_quals for v in allele_quals]
        else:
            weighted_base_freqs = [0] * len(ref_plus_alts)

        ref_qual = scaffold.ref_qual(pos) if ref else 0
        ref_fraction = round(scaffold.ref_fraction(pos), 3) if ref else 0.0

        info_dict = {
            'DP': int(scaffold.depth(pos)),
            'RQ': int(ref_qual),
            'MQ': int(round(scaffold.mean_mq(pos))),
            'RF': ref_fraction,
            'AD': allele_counts,
            'QS': allele_quals,
            'BF': weighted_base_freqs,
            'ST': [int(b in strong) for b in ref_plus_alts]
        }
        record = writer.new_record(
            scaffold.name,
            start=pos,
            alleles=[str(a) for a in ref_plus_alts],
            filter="PASS",
            info=info_dict,
            samples=[{"GT": "./."}],
        )

        yield record


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

    contig_lengths = []
    for scaffold in call_data.scaffolds_data.values():
        contig_lengths.append(
            f"##contig=<ID={scaffold.name},length={scaffold.length}>"
        )

    header = pysam.VariantHeader()
    header_str = VCF_TEMPLATE.format(
        date=datetime.now(),
        ref=call_data.reference_fasta,
        contig_lengths="\n".join(contig_lengths)
    )

    for line in header_str.split('\n'):
        header.add_line(line)

    header.add_sample("straingr")

    vcf_writer = pysam.VariantFile(output_file, 'w', header=header)

    record_iter = itertools.chain.from_iterable(
        vcf_records_for_scaffold(vcf_writer, scaffold, verboseness)
        for scaffold in call_data.scaffolds_data.values()
    )

    for record in record_iter:
        vcf_writer.write(record)
