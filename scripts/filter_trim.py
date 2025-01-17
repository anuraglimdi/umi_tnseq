#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to filter and retain reads containing transposon sequence at the 5' end.
Optionally trims reads to keep a specified number of bp after the transposon sequence

Input:
- fastq file with transposon sequencing data

Outputs:
- filtered fastq file (transpson sequence is clipped)
- list of unique molecular identifiers
- list of indices of reads that pass the filter 

Usage:
    python fastq_filter.py -i input.fastq -o output_directory [--trim_length N]

Example:
    python fastq_filter.py -i sample.fastq -o filtered_results --trim_length 50
"""

import time
import os
import argparse
import logging
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import numpy as np
import regex as re  # Using the 'regex' module for fuzzy matching


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Filter and optionally trim reads from a FastQ file based on a mariner sequence. "
        "Extracts UMIs and outputs indices of passing reads."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to the input fastq file (Note this should be the forward read (R1) that contains transposon sequence).",
    )
    parser.add_argument(
        "-o", "--output_dir", required=True, help="Directory to store the output files."
    )
    parser.add_argument(
        "--trim_length",
        type=int,
        default=None,
        help="Number of bases to retain after the mariner sequence match in the output FastQ.",
    )
    parser.add_argument(
        "--umi_length",
        type=int,
        default=10,
        help="Length of the UMI at the beginning of each read (default: 10).",
    )
    parser.add_argument(
        "--transposon_seq",
        type=str,
        default="GGGGACTTATCAGCCAACCTGTTA",
        help="Transposon sequence to search for (default: GGGGACTTATCAGCCAACCTGTTA).",
    )
    parser.add_argument(
        "--max_errors",
        type=int,
        default=1,
        help="Maximum number of allowed errors in the transposon sequence (default: 1).",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="filter_trim 1.0",
        help="Show the script version.",
    )
    return parser.parse_args()


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
    )


def ensure_output_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)
        logging.info(f"Created output directory: {path}")
    else:
        logging.info(f"Output directory already exists: {path}")


def filter_and_process_fastq(
    input_fastq, output_dir, transposon_seq, max_errors, trim_length, umi_length
):
    start_time = time.time()

    # Prepare output file paths
    base_filename = os.path.basename(input_fastq)
    trimmed_fastq = os.path.join(output_dir, f"filtered_{base_filename}")
    umi_file = os.path.join(output_dir, f"UMI_{os.path.splitext(base_filename)[0]}.txt")
    index_file = os.path.join(
        output_dir, f"PF_{os.path.splitext(base_filename)[0]}.index"
    )

    # Compile the regex pattern once for efficiency
    pattern = re.compile(f"({transposon_seq}){{e<={max_errors}}}", re.IGNORECASE)

    passed_indices = []
    counter = 0  # Overall read counter

    logging.info(f"Processing FastQ file: {input_fastq}")

    with open(input_fastq, "r") as in_handle, open(
        trimmed_fastq, "w"
    ) as out_handle, open(umi_file, "w") as umi_handle:

        for title, seq, qual in FastqGeneralIterator(in_handle):
            counter += 1
            match = pattern.search(seq)
            if match:
                passed_indices.append(counter - 1)  # Zero-based indexing

                # Determine the trim position
                trim_position = match.span()[1] - 2  # Retain 'TA' from original script

                if trim_length is not None:
                    # Calculate the end position ensuring it doesn't exceed the sequence length
                    end_position = trim_position + trim_length
                    trimmed_seq = seq[trim_position:end_position]
                    trimmed_qual = qual[trim_position:end_position]
                else:
                    # Retain the entire sequence after the trim position
                    trimmed_seq = seq[trim_position:]
                    trimmed_qual = qual[trim_position:]

                # Write the trimmed or original read to the output FastQ
                out_handle.write(f"@{title}\n{trimmed_seq}\n+\n{trimmed_qual}\n")

                # Extract and write the UMI
                umi = seq[:umi_length]
                umi_handle.write(f"{umi}\n")

            if counter % 1_000_000 == 0:
                logging.info(f"Processed {counter} reads...")

    # Save the indices of passing reads
    np.savetxt(index_file, passed_indices, fmt="%d")

    elapsed_time = time.time() - start_time
    logging.info(f"Filtering completed in {elapsed_time:.2f} seconds.")
    logging.info(f"Total reads processed: {counter}")
    logging.info(f"Total reads passed: {len(passed_indices)}")
    logging.info(f"Trimmed FastQ saved to: {trimmed_fastq}")
    logging.info(f"UMIs saved to: {umi_file}")
    logging.info(f"Indices saved to: {index_file}")


def main():
    setup_logging()
    args = parse_arguments()
    ensure_output_directory(args.output_dir)
    filter_and_process_fastq(
        input_fastq=args.input,
        output_dir=args.output_dir,
        transposon_seq=args.transposon_seq,
        max_errors=args.max_errors,
        trim_length=args.trim_length,
        umi_length=args.umi_length,
    )


if __name__ == "__main__":
    main()
