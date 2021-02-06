#!/usr/bin/env python3.8

import pysam
import argparse

import os
import sys
import time
import logging
import re
import statistics
import math

import subprocess
import tempfile

from collections import OrderedDict
from collections import namedtuple

################################################################################

# Min PL for an alignment to be kept as "good":
# This should be equal to the base read quality.
MIN_GOOD_ALIGNMENT_PL = 7.0

# Default max PL:
MAX_REPORTED_ALIGNMENT_PL = 60

# Minimum reported PL quality value (will override lower values):
MIN_REPORTED_ALIGNMENT_PL = 0

################################################################################


class ProcessedAlignmentResult(namedtuple(
    "ProcessedAlignmentResult", ["seq_name", "alignment_string", "target_start_index", "target_end_index",
                                 "read_start_pos", "read_end_pos", "template_length", "cigar", "overall_quality"]
)):
    __slots__ = ()

    def __str__(self):

        return (
            f"ProcessedAlignmentResult({self.seq_name}, {self.alignment_string}, {self.target_start_index}, "
            f"{self.target_end_index}, {self.read_start_pos}, {self.read_end_pos}, "
            f"{self.template_length}, {cigar_tuple_to_string(self.cigar)}, {self.overall_quality})"
        )

################################################################################


LOGGER = logging.getLogger("get_approx_raw_subread_array_lengths")


def configure_logging(args):
    """Set up logging for the module"""

    format_string = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'

    # Set logging level:
    log_level = logging.INFO
    if args.quiet:
        log_level = logging.CRITICAL
    elif args.verbose:
        log_level = logging.DEBUG
    elif args.veryverbose:
        log_level = logging.NOTSET

    logging.basicConfig(level=log_level, format=format_string)


################################################################################


def get_qual_pl(num_errors, seq_length):
    """
    Computes the Phred-Scaled quality score of an alignment.
    :param num_errors: Number of errors in the alignment.
    :param seq_length: Length of the alignment.
    :return: The PL quality.
    """

    if num_errors == 0:
        return MAX_REPORTED_ALIGNMENT_PL
    else:
        q = -10 * math.log10(num_errors / seq_length)
        if q < MIN_REPORTED_ALIGNMENT_PL:
            return MIN_REPORTED_ALIGNMENT_PL
        else:
            return q


def ingest_fastx_file(file_path):
    """Ingest the contents of a FASTA/FASTQ file and return a dictionary containing the sequence info"""

    _read_to_sequence_dict = OrderedDict()
    with pysam.FastxFile(file_path) as file_handle:
        for entry in file_handle:
            _read_to_sequence_dict[entry.name] = entry.sequence

    return _read_to_sequence_dict


def create_alignment_with_bwa_aln(read_name, read_sequence, delimiter_dict, minqual, threads=1, log_spacing="    "):
    """
    Perform a BWA ALN alignment on

    :param read_name: name of the read to which to align delimiters
    :param read_sequence: raw bases of a sequence to which to align the given delimiters
    :param delimiter_dict: Dictionary of delimiter sequences to align to the given sequence
    :param minqual: Minimum quality for an alignment to be retained.
    :param minbases: Minimum number of bases for an alignment to be retained.
    :param threads: number of threads to use when aligning.
    :param log_spacing: Spacing to precede any log statements.
    :return: A list of A list of ProcessedAlignmentResult objects. objects.
    """

    out_sai_file = "tmp.sai"
    out_file_name = "tmp.sam"

    # Write sequences to tmp fasta file:
    _, ref_file = tempfile.mkstemp()
    _, seq_file = tempfile.mkstemp()

    try:
        LOGGER.debug("%sCreating tmp \"reference\" fasta file: %s", log_spacing, ref_file)
        with open(ref_file, "w", encoding="ascii") as tmp:
            tmp.write(f">{read_name}\n")
            tmp.write(f"{read_sequence}\n")

        bwa_index_args = ["/bwa/bwa", "index", ref_file]
        LOGGER.debug(
            "%sRunning BWA Index on tmp file (%s): %s", log_spacing, ref_file, " ".join(bwa_index_args)
        )
        _ = subprocess.run(bwa_index_args, capture_output=True, check=True)

        LOGGER.debug("%sCreating tmp known segment \"read\" fasta file: %s", log_spacing, seq_file)
        seen_targets = set()
        with open(seq_file, "w", encoding="ascii") as tmp:
            for k,v in delimiter_dict.items():
                if k not in seen_targets:
                    tmp.write(f">{k}\n")
                    tmp.write(f"{v}\n")
                    seen_targets.add(k)

        LOGGER.debug("Contents of tmp \"reference\" fasta file:")
        with open(ref_file, "r") as f:
            for l in f.readlines():
                LOGGER.debug(l.rstrip())

        LOGGER.debug("Contents of known segment \"read\" fasta file:")
        with open(seq_file, "r") as f:
            for l in f.readlines():
                LOGGER.debug(l.rstrip())

        # Run the alignment:
        bwa_aln_args = [
            "/bwa/bwa", "aln",
            "-n 0.04",               # max #diff (int) or missing prob under 0.02 err rate (float) [0.04]
            "-o 1",                  # maximum number or fraction of gap opens [1]
            "-e -1",                 # maximum number of gap extensions, -1 for disabling long gaps [-1]
            "-d 10",                 # maximum occurrences for extending a long deletion [10]
            "-k 2",                  # maximum differences in the seed [2]
            "-M 3",                  # mismatch penalty [3]
            "-O 11",                 # gap open penalty [11]
            "-E 4",                  # gap extension penalty [4]
            "-l 8",                  # seed length [32]
            "-i 1",                  # do not put an indel within INT bp towards the ends [5]
            "-R 30",                 # stop searching when there are >INT equally best hits
            "-q 0",                  # quality threshold for read trimming down to 35bp [0]
            "-t", str(threads),
            "-f", out_sai_file,
            ref_file,
            seq_file
        ]

        LOGGER.debug(
            "%sRunning BWA ALN on tmp file (%s): %s", log_spacing, ref_file, " ".join(bwa_aln_args)
        )
        subprocess.run(bwa_aln_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True, text=True)

        # Convert the alignment to sam file:
        bwa_samse_args = [
            "/bwa/bwa", "samse",
            f"-n {len(target_sequences)}",  # Maximum number of alignments to output in the XA tag for
                                            # reads paired properly. If a read has more than INT hits, the
                                            # XA tag will not be written. [3]
            f"-f{out_file_name}",
            ref_file,
            out_sai_file,
            seq_file
        ]

        LOGGER.debug(
            f"{log_spacing}Running BWA SAMSE on tmp SAI file ({out_sai_file}): {' '.join(bwa_samse_args)}"
        )
        completed_process = subprocess.run(bwa_samse_args,
                                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True, text=True)

        if LOGGER.isEnabledFor(logging.DEBUG):
            LOGGER.debug("BWA ALN output:")
            for l in completed_process.stdout.split("\n"):
                LOGGER.debug(l)

            with open(out_file_name, 'rb') as f:
                LOGGER.debug("=" * 80)
                LOGGER.debug("Raw BWA ALN Alignment:")
                for line in f.readlines():
                    LOGGER.debug("%s", line.decode("ascii").rstrip())
                LOGGER.debug("=" * 80)

        return get_processed_results_from_bwa_aln_file(out_file_name, minqual)

    except subprocess.CalledProcessError as e:
        LOGGER.error("Could not align with BWA ALN!")
        LOGGER.error("Stdout: %s", e.stdout.decode("utf-8"))
        LOGGER.error("Stderr: %s", e.stderr.decode("utf-8"))
        raise e

    finally:
        os.remove(ref_file)
        os.remove(seq_file)
        try:
            os.remove(out_sai_file)
            pass
        except FileNotFoundError:
            # If the alignment failed, we won't necessarily have an output file:
            pass
        try:
            os.remove(out_file_name)
            pass
        except FileNotFoundError:
            # If the alignment failed, we won't necessarily have an output file:
            pass


def get_processed_results_from_bwa_aln_file(file_path, minqual):
    """
    Ingests the given file and creates raw alignments from it.

    :param file_path: Path to the sam file of a BWA ALN run.
    :param minqual: Minimum quality for an alignment to be retained.
    :return: A list of ProcessedAlignmentResult objects.
    """

    processed_results = []

    # Only primary hits will have sequence information with them so we have to
    # Read it off first.  There shouldn't be very many reads, so this is a little slow, but should be OK.
    read_seqs = dict()
    with pysam.AlignmentFile(file_path, 'r', check_sq=False) as f:
        for read in f.fetch(until_eof=True):
            if read.query_sequence:
                read_seqs[read.query_name] = read.query_sequence

    with pysam.AlignmentFile(file_path, 'r', check_sq=False) as f:
        for read in f.fetch(until_eof=True):

            if read.is_unmapped:
                continue

            seq_name = read.query_name
            bases = read_seqs[read.query_name]

            template_length = read.infer_read_length()
            qual_pl = get_qual_pl(read.get_tag("NM"), template_length)

            # Get the leading and trailing clips so we can remove them from the aligned string:
            leading_clips = 0
            for e in read.cigartuples:
                if e[0] == pysam.CSOFT_CLIP or e[0] == pysam.CHARD_CLIP:
                    leading_clips += e[1]
                else:
                    break

            trailing_clips = 0
            for e in read.cigartuples[::-1]:
                if e[0] == pysam.CSOFT_CLIP or e[0] == pysam.CHARD_CLIP:
                    trailing_clips += e[1]
                else:
                    break

            # Make a list of candidate read alignments so we can filter them all later:
            candidate_alignments = []

            # Note - must adjust ref start/end pos to align properly with conventions from other aligners
            #        (other aligner conventions: 1-based coordinates, inclusive end positions)
            p = ProcessedAlignmentResult(
                seq_name, bases, leading_clips, len(bases)-trailing_clips-1,
                int(read.reference_start), int(read.reference_end-1),
                template_length, tuple(read.cigartuples), qual_pl
            )
            candidate_alignments.append(p)

            # # Get all XA alignments as well:
            # candidate_alignments.extend(process_xa_tags(read, read_seqs))

            # Process all alignments with the same rules:
            for p in candidate_alignments:
                # Check against thresholds to make sure we should report the alignment:
                if p.overall_quality < minqual:
                    reason_string = f"qual too low ({p.overall_quality} < {minqual})"
                    LOGGER.debug("Target does not pass threshold: %s: %s (%s)", p.seq_name, reason_string, p)
                else:
                    processed_results.append(p)

    # Sort by the order in which they appear in the read.
    # This is _VERY_IMPORTANT_ for finding the ordered regions in a following step.
    processed_results.sort(key=lambda x: x.read_start_pos)

    for r in processed_results:
        LOGGER.debug("    %s", r)

    return processed_results


def get_array_element_alignments(args):
    """Get the stats from the given polymerase reads and dump them to the specified output file(s)."""

    # read in our sequence delimiters:
    mas_seq_delimiters = ingest_fastx_file(args.delimiters)

    stat_names = ['num_array_elements']

    # Print header:
    with open(args.outfile, 'w') as out_file:

        out_file.write("\t".join(['zmw'] + stat_names + ['array_composition', 'concise_alignment_info']))
        out_file.write("\n")

        with pysam.AlignmentFile(args.bam, 'rb', check_sq=False) as bam_file:

            zmw_re = re.compile(r'.*?/(.*?)/.*')

            prev_zmw = None
            zmw_read_lengths = []
            subread_seqs = []
            subread_names = []

            num_reads = 0

            for r in bam_file.fetch(until_eof=True):

                zmw = int(zmw_re.match(r.query_name).group(1))
                if prev_zmw and zmw != prev_zmw:

                    stats = dict()
                    stats["mean"] = statistics.mean(zmw_read_lengths)
                    stats["median"] = statistics.median(zmw_read_lengths)

                    read_length_diffs = list(
                        map(
                            lambda x: (abs(x - stats["mean"]) + abs(x - stats["median"]))/2,
                            zmw_read_lengths
                        )
                    )
                    index_min = min(range(len(read_length_diffs)), key=read_length_diffs.__getitem__)

                    ###########################################
                    # Analyze read here:
                    #######################
                    LOGGER.info(f"ZMW {zmw} \"Best\" subread: {subread_names[index_min]} -> {subread_seqs[index_min]}")

                    processed_results = create_alignment_with_bwa_aln(
                        r.query_name,
                        r.query_sequence,
                        mas_seq_delimiters,
                        minqual=args.minqual
                    )

                    #     "ProcessedAlignmentResult", ["seq_name", "alignment_string", "target_start_index", "target_end_index",
                    #                                  "read_start_pos", "read_end_pos", "template_length", "cigar", "overall_quality"]

                    concise_alignments = ",".join([
                            f"{p.seq_name}"
                            f"[@{p.read_start_pos}"
                            f"q{p.overall_quality}"
                            f"c{100.0*(p.target_end_index - p.target_start_index)/p.template_length:2.0f}]"
                            for p in processed_results
                        ])

                    # Write our output to the file:
                    out_file.write(
                        f"{zmw}\t"
                        f"{len(processed_results)}\t"
                        f"{','.join([p.seq_name for p in processed_results])}"
                        f"{concise_alignments}\n"
                    )

                    ###########################################
                    # Reset for next ZMW:
                    #######################
                    zmw_read_lengths = []
                    subread_seqs = []
                    subread_names = []

                zmw_read_lengths.append(len(r.query_sequence))
                subread_seqs.append(r.query_sequence)
                subread_names.append(r.query_name)
                prev_zmw = zmw

                num_reads += 1

                if num_reads % 5000 == 0:
                    LOGGER.info(f"Processed {num_reads} reads.")

            # Now we have to handle the last ZMW:
            # TODO: DO THIS

        LOGGER.info(f"Total reads processed: {num_reads}")

################################################################################


def main(raw_args):

    base_outfile_name = "approx_raw_subread_array_lengths"

    # Get our start time:
    overall_start = time.time()

    parser = argparse.ArgumentParser(
        description="Ingests two files: "
                    "A PacBio subreads file (SAM/BAM) containing raw subreads off the instrument."
                    "(Under the hood uses BWA-MEM 2)"
                    "A FASTA file containing the sequences of the delimiters in the MASSeq library",
        usage="Get approximate MASSeq library composition from raw PacBio subreads.",
    )

    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument(
        "-b", "--bam", help="PacBio subreads SAM/BAM file.", required=True
    )

    required_args.add_argument(
        "-d", "--delimiters", help="FASTA file containing MASSeq delimiters.", required=True
    )

    parser.add_argument(
        "-m",
        "--minqual",
        help=f"Minimum quality for good alignment (default: {MIN_GOOD_ALIGNMENT_PL}).  "
             f"Generally this should be the same as the average base quality.",
        default=MIN_GOOD_ALIGNMENT_PL,
        type=float
    )

    parser.add_argument(
        "-o",
        "--outfile",
        help="Output file in which to store results.",
        default=f"{base_outfile_name}.tsv",
        required=False,
    )

    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument(
        "-q", "--quiet", help="silence logging except errors", action="store_true"
    )
    verbosity_group.add_argument(
        "-v", "--verbose", help="increase output verbosity", action="store_true"
    )
    verbosity_group.add_argument(
        "-vv", "--veryverbose", help="maximal output verbosity", action="store_true"
    )

    # ---------------------------------------

    # Parse args
    args = parser.parse_args()

    configure_logging(args)

    # Log our command-line and log level so we can have it in the log file:
    LOGGER.info("Invoked by: %s", " ".join(raw_args))
    LOGGER.info("Complete runtime configuration settings:")
    for name, val in vars(args).items():
        LOGGER.info("    %s = %s", name, val)
    LOGGER.info("Log level set to: %s", logging.getLevelName(logging.getLogger().level))

    # Call our main method:
    get_array_element_alignments(args)

    LOGGER.info(f"Results written to: {args.outfile}")

    overall_end = time.time()
    LOGGER.info("Elapsed time: %f", overall_end - overall_start)


################################################################################


if __name__ == '__main__':
    main(sys.argv)
