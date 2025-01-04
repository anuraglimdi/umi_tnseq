#!/usr/bin/env bash
# =============================================================================
# Script: run_bowtie.sh
# Description: runs bowtie2, an alignment tool for short reads on a single FastQ file, 
#              aligns reads to a reference genome and outputs a SAM file to a specified directory.
# =============================================================================

set -euo pipefail  # Enable strict error handling

# -----------------------------------
# Function: print_usage
# Description: Display script usage instructions
# -----------------------------------
print_usage() {
    echo "Usage: $(basename "$0") -i INPUT_FASTQ -o OUTPUT_DIR -r REFERENCE_GENOME [-p THREADS]"
    echo ""
    echo "Options:"
    echo "  -i, --input        Path to the input FastQ file."
    echo "  -o, --output-dir   Directory to store the output SAM file."
    echo "  -r, --reference    Prefix of the Bowtie2 index for the reference genome."
    echo "  -p, --threads      (Optional) Number of threads to use (default: 1)."
    echo "  -h, --help         Display this help message and exit."
    echo ""
    echo "Example:"
    echo "  $(basename "$0") -i sample.fastq -o aligned_reads -r genome_index -p 4"
}

# -----------------------------------
# Function: log
# Description: Log messages with timestamp
# -----------------------------------
log() {
    local type="$1"
    local message="$2"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$type] $message"
}

# -----------------------------------
# Parse command-line arguments
# -----------------------------------
# Initialize default values
THREADS=1

# Use getopt for parsing long and short options
PARSED_ARGS=$(getopt -o i:o:r:p:h --long input:,output-dir:,reference:,threads:,help -- "$@")
if [[ $? -ne 0 ]]; then
    print_usage
    exit 1
fi

eval set -- "$PARSED_ARGS"

# Extract options and their arguments into variables
while true; do
    case "$1" in
        -i|--input)
            INPUT_FASTQ="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -r|--reference)
            REFERENCE_GENOME="$2"
            shift 2
            ;;
        -p|--threads)
            THREADS="$2"
            shift 2
            ;;
        -h|--help)
            print_usage
            exit 0
            ;;
        --)
            shift
            break
            ;;
        *)
            log "ERROR" "Unknown option: $1"
            print_usage
            exit 1
            ;;
    esac
done

# Check mandatory arguments
if [[ -z "${INPUT_FASTQ:-}" ]] || [[ -z "${OUTPUT_DIR:-}" ]] || [[ -z "${REFERENCE_GENOME:-}" ]]; then
    log "ERROR" "Missing required arguments."
    print_usage
    exit 1
fi

# -----------------------------------
# Validate input FastQ file
# -----------------------------------
if [[ ! -f "$INPUT_FASTQ" ]]; then
    log "ERROR" "Input FastQ file does not exist: $INPUT_FASTQ"
    exit 1
fi

# -----------------------------------
# Validate Bowtie2 reference index
# -----------------------------------
# Bowtie2 expects index files with extensions like .1.bt2, .2.bt2, etc.
# We'll check for the existence of the first index file as a proxy
if [[ ! -f "${REFERENCE_GENOME}.1.bt2" ]]; then
    log "ERROR" "Bowtie2 reference index not found: ${REFERENCE_GENOME}.1.bt2"
    log "INFO" "Ensure that the Bowtie2 index is correctly built with prefix '${REFERENCE_GENOME}'."
    exit 1
fi

# -----------------------------------
# Prepare output directory
# -----------------------------------
if [[ ! -d "$OUTPUT_DIR" ]]; then
    mkdir -p "$OUTPUT_DIR"
    log "INFO" "Created output directory: $OUTPUT_DIR"
else
    log "INFO" "Output directory exists: $OUTPUT_DIR"
fi

# -----------------------------------
# Define output SAM file path
# -----------------------------------
# Extract base name without path and extension
BASE_NAME=$(basename "$INPUT_FASTQ")
BASE_NAME="${BASE_NAME%%.*}"  # Removes extension (first dot onwards)

OUTPUT_SAM="${OUTPUT_DIR}/${BASE_NAME}.sam"

# -----------------------------------
# Run Bowtie2
# -----------------------------------
log "INFO" "Starting alignment:"
log "INFO" "Input FastQ: $INPUT_FASTQ"
log "INFO" "Reference Genome Index: $REFERENCE_GENOME"
log "INFO" "Output SAM: $OUTPUT_SAM"
log "INFO" "Threads: $THREADS"

# Execute Bowtie2
bowtie2 -q \
        -x "$REFERENCE_GENOME" \
        -U "$INPUT_FASTQ" \
        -S "$OUTPUT_SAM" \
        -p "$THREADS" \
        --silent 2> "${OUTPUT_DIR}/${BASE_NAME}_bowtie2.log"

log "INFO" "Bowtie2 alignment completed successfully."

# Optionally, you can convert SAM to BAM, sort, and index using samtools
# Uncomment the following lines if desired

# log "INFO" "Converting SAM to BAM..."
# samtools view -bS "$OUTPUT_SAM" > "${OUTPUT_SAM%.sam}.bam"
# log "INFO" "Sorting BAM file..."
# samtools sort "${OUTPUT_SAM%.sam}.bam" -o "${OUTPUT_SAM%.sam}_sorted.bam"
# log "INFO" "Indexing BAM file..."
# samtools index "${OUTPUT_SAM%.sam}_sorted.bam"
# log "INFO" "BAM processing completed."

exit 0