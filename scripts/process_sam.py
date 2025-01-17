#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: process_sam.py
Description: Processes paired SAM files (Read1 and Read2), integrates UMI data,
             discards PCR duplicates, and outputs the number of reads for each TA site.
Usage:
    python process_sam.py \
        --sam_r1 /path/to/read1.sam \
        --sam_r2 /path/to/read2.sam \
        --umi_list /path/to/UMI_sample.txt \
        --indices_pf /path/to/PF_sample.index \
        --output_dir /path/to/output_directory
"""

import os
import sys
import argparse
import pysam
import numpy as np
import logging
import re


def parse_arguments():
    """
    Parses command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Process paired SAM files (Read1 and Read2), integrate UMI data, \
                     discard PCR duplicates, and output read counts per TA site."
    )
    parser.add_argument("--sam_r1", required=True, help="Path to the Read1 SAM file.")
    parser.add_argument("--sam_r2", required=True, help="Path to the Read2 SAM file.")
    parser.add_argument("--umi_list", required=True, help="Path to the UMI list file.")
    parser.add_argument(
        "--indices_pf", required=True, help="Path to the PF index file."
    )
    parser.add_argument(
        "--output_dir", required=True, help="Directory to store the output file."
    )
    return parser.parse_args()


def setup_logging(output_dir):
    """
    Sets up logging configuration.
    """
    log_file = os.path.join(output_dir, "process_sam.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
    )
    logging.info("Logging initialized.")


def get_positions(sam_file_path):
    """
    Extracts mapping positions from the SAM file.
    Returns a NumPy array of positions.
    """
    positions = []
    try:
        with pysam.AlignmentFile(sam_file_path, "r") as sam_file:
            for read in sam_file:
                pos = read.get_reference_positions(full_length=False)
                # Check if the read is mapped and uniquely aligned
                if pos and not read.has_tag("XS"):
                    # Determine strand
                    if read.is_forward:
                        positions.append(pos[0])
                    elif read.is_reverse:
                        positions.append(pos[-1])
                    else:
                        positions.append(-1)  # Undefined strand
                else:
                    positions.append(-1)  # Unmapped or multiple mappings
    except Exception as e:
        logging.error(f"Error processing SAM file {sam_file_path}: {e}")
        sys.exit(1)

    positions = np.array(positions)
    logging.info(f"Extracted positions from SAM file: {sam_file_path}")
    return positions


def load_umi_list(umi_file):
    """
    Loads UMI list from a text file.
    Returns a list of UMIs.
    """
    if not os.path.isfile(umi_file):
        logging.error(f"UMI file does not exist: {umi_file}")
        sys.exit(1)
    try:
        with open(umi_file, "r") as f:
            umi_list = f.read().splitlines()
        logging.info(f"Loaded {len(umi_list)} UMIs from {umi_file}")
        return umi_list
    except Exception as e:
        logging.error(f"Error reading UMI file {umi_file}: {e}")
        sys.exit(1)


def load_indices_pf(indices_file):
    """
    Loads indices of reads passing filter from a file.
    Returns a NumPy array of indices.
    """
    if not os.path.isfile(indices_file):
        logging.error(f"Indices PF file does not exist: {indices_file}")
        sys.exit(1)
    try:
        indices = np.loadtxt(indices_file, dtype=int)
        logging.info(f"Loaded {len(indices)} PF indices from {indices_file}")
        return indices
    except Exception as e:
        logging.error(f"Error loading indices PF file {indices_file}: {e}")
        sys.exit(1)


def process_umi_positions(umi_list, positions_r2_pf):
    """
    Combines UMI with Read2 mapping positions to create unique identifiers.
    Returns a list of combined UMI+position strings.
    """
    if len(umi_list) != len(positions_r2_pf):
        logging.error(
            "Length of UMI list does not match length of Read2 positions after filtering."
        )
        sys.exit(1)
    combined_umi = [f"{umi}{pos}" for umi, pos in zip(umi_list, positions_r2_pf)]
    logging.info("Combined UMI with Read2 positions.")
    return combined_umi


def discard_pcr_duplicates(positions_combined, combined_umi):
    """
    Discards PCR duplicates based on unique UMIs per position.
    Returns arrays of unique positions, raw counts, and UMI-corrected counts.
    """
    # Ensure non-negative positions
    valid_indices = positions_combined != -1
    filtered_positions = positions_combined[valid_indices]
    filtered_umi = [
        umi for umi, pos in zip(combined_umi, positions_combined) if pos != -1
    ]

    # Unique positions and counts
    unique, counts = np.unique(filtered_positions, return_counts=True)

    # Initialize UMI-corrected counts
    counts_umi = np.zeros(len(unique), dtype=int)

    # Sort UMIs by positions to facilitate grouping
    sorted_indices = np.argsort(filtered_positions)
    sorted_positions = filtered_positions[sorted_indices]
    sorted_umis = [filtered_umi[i] for i in sorted_indices]

    # Group UMIs by position and count unique UMIs
    unique_sorted, cumulative = np.unique(sorted_positions, return_counts=True)
    cumulative = np.cumsum(cumulative)
    start = 0
    for i in range(len(unique_sorted)):
        end = cumulative[i]
        umi_set = set(sorted_umis[start:end])
        counts_umi[i] = len(umi_set)
        start = end

    logging.info("Discarded PCR duplicates based on UMI data.")
    return unique_sorted, counts, counts_umi


def count_reads_per_ta_site(unique_positions, counts, counts_umi):
    """
    Constructs a matrix with:
    - Column 0: Unique Positions
    - Column 1: Raw Counts
    - Column 2: UMI-corrected Counts
    Returns the matrix as a NumPy array.
    """
    position_counts = np.vstack((unique_positions, counts, counts_umi)).astype(int)
    logging.info("Constructed position counts matrix.")
    return position_counts


def save_position_counts(matrix, output_dir, sample_name):
    """
    Saves the position counts matrix to a file.
    """
    output_file = os.path.join(output_dir, f"{sample_name}_merged.pos")
    try:
        header = "Position\tRaw_Count\tUMI_Corrected_Count"
        np.savetxt(
            output_file, matrix.T, fmt="%d", delimiter="\t", header=header, comments=""
        )
        logging.info(f"Saved position counts to {output_file}")
    except Exception as e:
        logging.error(f"Failed to save position counts to {output_file}: {e}")
        sys.exit(1)


def extract_sample_name(sam_r1_path):
    """
    Extracts the sample name from the read 1 SAM filename.
    Assumes the filename is in the format: {sample}_R1.sam
    Adjust the regex pattern if your filenames follow a different convention.
    """
    sam_r1_filename = os.path.basename(sam_r1_path)
    match = re.match(r".*/(?P<sample_name>.+)_R1_.*$", sam_r1_filename)
    if match:
        return match.group("sample")
    else:
        logging.error(
            f"Read 1 SAM filename does not match expected pattern: {sam_r1_filename}"
        )
        sys.exit(1)


def main():
    # Parse command-line arguments
    args = parse_arguments()

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # Setup logging
    setup_logging(args.output_dir)

    # Extract sample name from Read1 SAM filename
    sample_name = extract_sample_name(args.sam_r1)
    logging.info(f"Processing sample: {sample_name}")

    # Extract positions from SAM files
    positions_r1 = get_positions(args.sam_r1)
    positions_r2 = get_positions(args.sam_r2)

    # Load UMI list
    umi_list = load_umi_list(args.umi_list)

    # Load PF indices
    indices_pf = load_indices_pf(args.indices_pf)

    # Slice out the data for Read2 corresponding to reads which passed the filter in Read1
    positions_r2_pf = positions_r2[indices_pf]
    umi_list_pf = [umi_list[x] for x in indices_pf]

    # Assert that the length of the Read1 and Read2 PF positions array is equal
    # Doing this check because we have already filtered out invalid reads without transposon
    # for read 1 in the first step of the analysis
    if len(positions_r1) != len(positions_r2_pf):
        logging.error(
            "Number of reads in read 1 and read 2 data is not equal after filtering."
        )
        sys.exit(1)

    # Create a list of UMIs where we incorporate the Read2 mapping coordinate
    combined_umi = process_umi_positions(umi_list_pf, positions_r2_pf)

    # Assert that the combined positions and the UMI list have the same length
    if len(positions_r1) != len(combined_umi):
        logging.error("Number of UMIs and mapped reads is not equal.")
        sys.exit(1)

    # Counting how many unique data points we have, and corresponding counts
    unique_positions, counts, counts_umi = discard_pcr_duplicates(
        positions_r1, combined_umi
    )

    # Construct the final matrix
    # Column 0 - list of all the unique positions
    # Column 1 - counts corresponding to those positions
    # Column 2 - UMI corrected counts corresponding to those positions
    position_counts_matrix = count_reads_per_ta_site(
        unique_positions, counts, counts_umi
    )

    # Save the position counts
    save_position_counts(position_counts_matrix, args.output_dir, sample_name)

    logging.info("Processing completed successfully.")


if __name__ == "__main__":
    main()
