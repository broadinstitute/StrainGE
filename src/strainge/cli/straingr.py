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
import sys
import logging
import argparse
from pathlib import Path

import pysam

from strainge.variant_caller import VariantCaller, Reference
from strainge.sample_compare import SampleComparison
from strainge.io.variants import (call_data_from_hdf5, call_data_to_hdf5,
                                  boolean_array_to_bedfile,  write_vcf,
                                  generate_call_summary_tsv, array_to_bedgraph)
from strainge.io.comparisons import (generate_compare_summary_tsv,
                                     generate_compare_details_tsv)
from strainge.cli.registry import Subcommand

logger = logging.getLogger()


def write_tracks(call_data, track_covered=None, track_coverage=None,
                 track_poor_mq=None, track_lowmq_count=None,
                 track_high_coverage=None, track_gaps=None, track_min_size=1):
    if track_covered:
        logger.info("Writing 'callable' BED track...")
        for scaffold in call_data.scaffolds_data.values():
            boolean_array_to_bedfile(scaffold.strong > 0, track_covered,
                                     scaffold.name, track_min_size)

    if track_coverage:
        logger.info("Writing 'coverage' BedGraph track...")

        # Write BedGraph trackline
        print("track type=bedGraph", file=track_coverage)

        # Write BedGraph values
        for scaffold in call_data.scaffolds_data.values():
            array_to_bedgraph(scaffold.coverage, track_coverage, scaffold.name)

    if track_lowmq_count:
        logger.info("Writing 'lowmq count' BedGraph track...")

        # Write BedGraph trackline
        print("track type=bedGraph", file=track_lowmq_count)

        # Write BedGraph values
        for scaffold in call_data.scaffolds_data.values():
            array_to_bedgraph(scaffold.lowmq_count, track_lowmq_count,
                              scaffold.name)

    if track_poor_mq:
        logger.info("Writing 'poor mapping quality' BED track...")
        for scaffold in call_data.scaffolds_data.values():
            boolean_array_to_bedfile(scaffold.lowmq, track_poor_mq,
                                     scaffold.name, track_min_size)

    if track_high_coverage:
        logger.info("Writing 'high coverage' BED track...")
        for scaffold in call_data.scaffolds_data.values():
            boolean_array_to_bedfile(scaffold.high_coverage,
                                     track_high_coverage,
                                     scaffold.name, track_min_size)

    if track_gaps:
        logger.info("Writing 'gaps' BED track...")
        writer = csv.writer(track_gaps, delimiter='\t', lineterminator='\n')
        for scaffold in call_data.scaffolds_data.values():
            for gap in scaffold.gaps:
                writer.writerow((scaffold.name, gap.start, gap.end))


class CallSubcommand(Subcommand):
    """
    StrainGR: strain-aware variant caller for metagenomic samples

    This command analyzes
    """
    def register_arguments(self, subparser: argparse.ArgumentParser):
        subparser.add_argument(
            'reference',
            help="Reference FASTA file. Can be GZIP compressed."
        )

        subparser.add_argument(
            'sample',
            help="BAM file with the aligned reads of the sample against the "
                 "reference"
        )

        call_qc_group = subparser.add_argument_group(
            'Quality control',
            'Options which determine which reads to consider, when base or '
            'read mapping qualities are high enough for calling, etc.'
        )
        call_qc_group.add_argument(
            '-Q', '--min-qual', type=int, default=5,
            help="Minimum quality for a base to be considered. Default: %("
                 "default)d."
        )
        call_qc_group.add_argument(
            '-P', '--min-pileup-qual', type=int, default=50,
            help="Minimum sum of qualities for an allele to be trusted."
                 "for variant calling. Default: %(default)d."
        )
        call_qc_group.add_argument(
            '-F', '--min-qual-frac', type=float, default=0.1,
            help="Minimum fraction of the reads in the pileup required to "
                 "confirm an allele (fractions are base quality weighted). "
                 "Default: %(default)g"
        )
        call_qc_group.add_argument(
            '-M', '--min-mapping-qual', type=int, default=5,
            help="Minimum mapping quality of the whole read to be considered. "
                 "Default: %(default)d."
        )
        call_qc_group.add_argument(
            '-G', '--min-gap', type=int, default=2000,
            help="Minimum size of gap to be considered as such. Default: "
                 "2000. Will be automatically scaled depending on coverage."
        )

        call_out_group = subparser.add_argument_group(
            "Output formats",
            "Options for writing the results to different file "
            "formats. If you're unsure what to choose, output at least the "
            "data to HDF5, the rest of the output files can be created "
            "afterwards from the HDF5 data using 'straingr view'."
        )

        call_out_group.add_argument(
            '-s', '--summary', type=argparse.FileType('w'), required=False,
            metavar='FILE', default=sys.stdout,
            help="Output a TSV with a summary of variant calling statistics "
                 "to the given file."
        )
        call_out_group.add_argument(
            '--hdf5-out', default=None, required=False, metavar='FILE',
            help="Output StrainGR variant calling data to the given HDF5 file."
        )
        call_out_group.add_argument(
            '-V', '--vcf', type=argparse.FileType('w'), required=False,
            default=None, metavar='FILE',
            help="Output a VCF file with SNP's. Please be aware that "
                 "we do not have a good insertion/deletion calling mechanism, "
                 "but some information on possible indels is written to the "
                 "VCF file."
        )
        call_out_group.add_argument(
            '--verbose-vcf', action="store_true", default=False,
            help="To be used with --vcf. If you set this flag, then the VCF "
                 "will also include records for positions in the genome where "
                 "nothing but the reference base is observed. By default it "
                 "will only output records for positions where some evidence "
                 "for a SNP is observed."
        )
        call_out_group.add_argument(
            '--track-callable', type=argparse.FileType('w'), required=False,
            default=None, metavar='FILE',
            help="Output a BED file to the given filename, indicating the "
                 "callable regions in the genome."
        )
        call_out_group.add_argument(
            '--track-coverage', type=argparse.FileType('w'), required=False,
            default=None, metavar='FILE',
            help="Output a BedGraph file to the given filename, containing "
                 "the coverage per position as seen by StrainGR."
        )
        call_out_group.add_argument(
            '--track-poor-mq', type=argparse.FileType('w'), required=False,
            default=None, metavar='FILE',
            help="Output a BED file to the given filename, indicating regions "
                 "with a majority of poorly mapped reads."
        )
        call_out_group.add_argument(
            '--track-lowmq-count', type=argparse.FileType('w'), required=False,
            default=None, metavar='FILE',
            help="Output a BedGraph file indicating how many reads with low "
                 "mapping quality are located at each position. This "
                 "includes counts from alternative alignment locations."
        )
        call_out_group.add_argument(
            '--track-high-coverage', type=argparse.FileType('w'),
            required=False,
            default=None, metavar='FILE',
            help="Output a BED file to the given filename, indicating regions "
                 "with abnormally high coverage."
        )
        call_out_group.add_argument(
            '--track-gaps', type=argparse.FileType('w'), required=False,
            default=None, metavar='FILE',
            help="Output a BED file to the given filename, indicating regions "
                 "marked as a gap"
        )
        call_out_group.add_argument(
            '--track-min-size', type=int, required=False, default=1,
            help="For all --track-* options above, only include features ("
                 "regions) of at least the given size. Default: %(default)d."
        )

    def __call__(self, reference, sample,
                 min_qual, min_pileup_qual, min_qual_frac,
                 min_mapping_qual, min_gap,
                 summary=None, hdf5_out=None,
                 vcf=None, verbose_vcf=False,
                 track_covered=None, track_coverage=None,
                 track_poor_mq=None, track_lowmq_count=None,
                 track_high_coverage=None,
                 track_gaps=None, track_min_size=1, **kwargs):
        """Call variants in a mixed-strain sample."""

        logger.info("Loading reference %s...", reference)
        reference = Reference(reference)
        logger.info("Reference length: %d", reference.length)
        sample_bam = pysam.AlignmentFile(sample)

        logger.info("Start analyzing aligned reads...")
        caller = VariantCaller(min_qual, min_pileup_qual, min_qual_frac,
                               min_mapping_qual, min_gap)

        call_data = caller.process(reference, sample_bam.pileup())

        if hdf5_out:
            # Output call datasets to HDF5
            logger.info("Writing data to HDF5 file %s...", hdf5_out)
            call_data_to_hdf5(call_data, hdf5_out)

        if summary:
            # Output a summary TSV
            if summary != sys.stdout:
                logger.info("Writing summary to %s", summary.name)

            generate_call_summary_tsv(call_data, summary)

        if vcf:
            logger.info("Generating VCF file...")
            write_vcf(call_data, vcf, not verbose_vcf)

        write_tracks(call_data, track_covered, track_coverage,
                     track_poor_mq, track_lowmq_count,
                     track_high_coverage, track_gaps, track_min_size)

        logger.info("Done.")


class ViewSubcommand(Subcommand):
    """
    View call statistics stored in a HDF5 file and output results to
    different file formats
    """

    def register_arguments(self, subparser: argparse.ArgumentParser):
        subparser.add_argument(
            'reference',
            help="Reference FASTA file. Can be GZIP compressed."
        )
        subparser.add_argument(
            'hdf5',
            help="HDF5 file with StrainGR call statistics."
        )

        subparser.add_argument(
            '-s', '--summary', type=argparse.FileType('w'), required=False,
            default=sys.stdout, metavar='FILE',
            help="Output a TSV with a summary of variant calling statistics "
                 "to the given file."
        )
        subparser.add_argument(
            '--track-covered', type=argparse.FileType('w'), required=False,
            default=None, metavar='FILE',
            help="Output a BED file to the given filename, indicating "
                 "callable regions in the genome."
        )
        subparser.add_argument(
            '--track-coverage', type=argparse.FileType('w'), required=False,
            default=None, metavar='FILE',
            help="Output a BedGraph file to the given filename, containing "
                 "the coverage per position as seen by StrainGR."
        )
        subparser.add_argument(
            '--track-poor-mq', type=argparse.FileType('w'), required=False,
            default=None, metavar='FILE',
            help="Output a  file to the given filename, indicating regions "
                 "with a majority of poorly mapped reads."
        )
        subparser.add_argument(
            '--track-lowmq-count', type=argparse.FileType('w'), required=False,
            default=None, metavar='FILE',
            help="Output a BedGraph file indicating how many reads with low "
                 "mapping quality are located at each position. This "
                 "includes counts from alternative alignment locations."
        )
        subparser.add_argument(
            '--track-high-coverage', type=argparse.FileType('w'),
            required=False,
            default=None, metavar='FILE',
            help="Output a BED file to the given filename, indicating regions "
                 "marked as abnormally high coverage."
        )
        subparser.add_argument(
            '--track-gaps', type=argparse.FileType('w'), required=False,
            default=None, metavar='FILE',
            help="Output a BED file to the given filename, indicating regions "
                 "marked as a gap"
        )
        subparser.add_argument(
            '--track-min-size', type=int, required=False, default=1,
            help="For all --track-* options above, only include features ("
                 "regions) of at least the given size. Default: %(default)d."
        )

        subparser.add_argument(
            '-V', '--vcf', type=argparse.FileType('w'), required=False,
            default=None, metavar='FILE',
            help="Output a VCF file with SNP's. Please be aware that "
                 "we do not have a good insertion/deletion calling mechanism, "
                 "but some information on possible indels is written to the "
                 "VCF file."
        )
        subparser.add_argument(
            '--verbose-vcf', action="store_true", default=False,
            help="To be used with --vcf. If you set this flag, then the VCF "
                 "will also include rows for positions in the genome where "
                 "nothing but the reference base is observed. By default it "
                 "will only output rows for positions where some evidence for "
                 "a SNP is observed."
        )

    def __call__(self, reference, hdf5, summary=None,
                 track_covered=None, track_coverage=None,
                 track_poor_mq=None, track_lowmq_count=None,
                 track_high_coverage=None, track_gaps=None, track_min_size=1,
                 vcf=None, verbose_vcf=False, **kwargs):
        """View and output the StrainGR calling results in different file
        formats."""
        logger.info("Loading reference %s...", reference)
        reference = Reference(reference)

        logger.info("Loading data from HDF5 file %s", hdf5)
        call_data = call_data_from_hdf5(reference, hdf5)

        if summary:
            # Output a summary TSV
            if summary != sys.stdout:
                logger.info("Writing summary to %s", summary.name)
            generate_call_summary_tsv(call_data, summary)

        write_tracks(call_data, track_covered, track_coverage,
                     track_poor_mq, track_lowmq_count,
                     track_high_coverage, track_gaps, track_min_size)

        if vcf:
            logger.info("Generating VCF file...")
            write_vcf(call_data, vcf, not verbose_vcf)

        logger.info("Done.")


class CompareSubCommand(Subcommand):
    """
    Compare strains and variant calls in two different samples. Reads of
    both samples must be aligned to the same reference.

    It's possible to generate a TSV with summary stats as well as a file
    with more detailed information on which alleles are called at what
    positions.
    """

    def register_arguments(self, subparser: argparse.ArgumentParser):
        subparser.add_argument(
            'reference', metavar='REF',
            help="The reference used for variant calling."
        )

        subparser.add_argument(
            'samples', nargs='+', metavar='SAMPLE_HDF5', type=Path,
            help="HDF5 files with variant calling data for each sample. "
                 "Number of samples should be exactly two, except when used "
                 "with --baseline."
        )

        subparser.add_argument(
            '-o', '--summary-out', type=argparse.FileType('w'),
            default=sys.stdout,
            help="Output file for summary statistics. Defaults to stdout."
        )

        subparser.add_argument(
            '-d', '--details-out', type=argparse.FileType('w'), default=None,
            help="Output file for detailed base level differences between "
                 "samples (optional)."
        )
        subparser.add_argument(
            '-V', '--verbose-details', action="store_true", default=False,
            help="Output detailed information for every position in the "
                 "genome instead of only for positions where alleles differ."
        )

        subparser.add_argument(
            '-b', '--baseline', default="", required=False, type=Path,
            help="Path to a sample to use as baseline, and compare all other "
                 "given samples to this one. Outputs a shell script that "
                 "runs all individual pairwise comparisons."
        )

        subparser.add_argument(
            '-D', '--output-dir', default="", required=False, type=Path,
            help="The output directory of all comparison files when using "
                 "--baseline."
        )

    def __call__(self, reference, samples, summary_out=None, details_out=None,
                 verbose_details=False, baseline=None, output_dir="",
                 *args, **kwargs):
        if baseline and not baseline.is_file() and not baseline == Path(""):
            logger.error("Baseline %s does not exists.", baseline)
            return 1
        elif baseline and baseline.is_file():
            output_dir = Path(output_dir)

            for sample in samples:
                if sample == baseline:
                    continue

                fname_base = f"{baseline.stem}.vs.{sample.stem}"
                summary_file = output_dir / f"{fname_base}.summary.tsv"
                details_file = output_dir / f"{fname_base}.details.tsv"
                print(sys.argv[0], "compare",
                      "-o", summary_file,
                      "-d", details_file,
                      reference, baseline, sample)
        else:
            if len(samples) != 2:
                logger.error("The number of samples given should be exactly "
                             "two. To compare multiple samples against a "
                             "single baseline, use --baseline.")

                return 1

            logger.info("Loading reference %s...", reference)
            reference = Reference(reference)

            logger.info("Comparing sample %s vs %s", samples[0].stem,
                        samples[1].stem)

            logger.info("Loading sample 1 %s", samples[0].stem)
            call_data1 = call_data_from_hdf5(reference, samples[0])
            logger.info("Loading sample 2 %s", samples[1].stem)
            call_data2 = call_data_from_hdf5(reference, samples[1])

            comparison = SampleComparison(call_data1, call_data2)

            generate_compare_summary_tsv(comparison, summary_out)

            if details_out:
                logger.info("Generating details file...")
                generate_compare_details_tsv(details_out, call_data1,
                                             call_data2, verbose_details)

            logger.info("Done.")
