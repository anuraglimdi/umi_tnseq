#!/bin/bash

# ------------------------------------------------------------
# Script Name: rename_fastq.sh
# Description: Renames FASTQ files from {sample}_R{1,2}_L{lane}.fastq
#              to {sample}_L{lane}_R{1,2}.fastq within a specified directory.
# Usage:       ./rename_fastq.sh /path/to/directory
# ------------------------------------------------------------

# Exit immediately if a command exits with a non-zero status
set -e

# Function to display usage information
usage() {
    echo "Usage: $0 /path/to/directory"
    exit 1
}

# Check if exactly one argument (directory path) is provided
if [ "$#" -ne 1 ]; then
    echo "Error: Exactly one directory path must be provided."
    usage
fi

# Assign the first argument to DIR variable
DIR="$1"

# Check if the provided argument is a valid directory
if [ ! -d "$DIR" ]; then
    echo "Error: '$DIR' is not a valid directory."
    usage
fi

# Iterate over all matching FASTQ files in the specified directory
shopt -s nullglob  # Enable nullglob to handle cases with no matching files
for file in "$DIR"/*_R[12]_L*.fastq; do
    # Extract the base filename (without the directory path)
    base_file=$(basename "$file")
    
    # Use regex to capture sample_name, read number, and lane identifier
    if [[ "$base_file" =~ ^(.+)_R([12])_L([0-9]{3})\.fastq$ ]]; then
        sample="${BASH_REMATCH[1]}"
        read_num="${BASH_REMATCH[2]}"
        lane="${BASH_REMATCH[3]}"
        
        # Construct the new filename with the desired format
        new_filename="${sample}_L${lane}_R${read_num}.fastq"
        
        # Full path for the new filename
        new_file="${DIR}/${new_filename}"
        
        # Perform the renaming operation
        mv "$file" "$new_file"
        
        # Output the renaming action
        echo "Renamed: '$base_file' -> '$new_filename'"
    else
        # Inform about files that do not match the expected pattern
        echo "Skipped: '$base_file' does not match the expected pattern."
    fi
done

# Inform the user if no files were found and renamed
shopt -u nullglob  # Disable nullglob
matched_files=("$DIR"/*_R[12]_L*.fastq)
if [ ${#matched_files[@]} -eq 0 ]; then
    echo "No FASTQ files matching the pattern were found in '$DIR'."
fi

exit 0