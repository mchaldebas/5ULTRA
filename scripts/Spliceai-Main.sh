#!/bin/bash

set -euo pipefail

# Function to display usage instructions
usage() {
    echo "Usage: $0 --input input_file [--output output_file]"
    exit 1
}

# Default values
input_file=""
output_file=""
intervals="./data/5UTRs.intervals.bed"

# Parse command-line arguments
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case $1 in
        --input)
            if [[ -n ${2-} ]]; then
                input_file="$2"
                shift 2
            else
                echo "Error: --input requires a non-empty argument."
                usage
            fi
            ;;
        --output)
            if [[ -n ${2-} ]]; then
                output_file="$2"
                shift 2
            else
                echo "Error: --output requires a non-empty argument."
                usage
            fi
            ;;
        -*|--*)
            echo "Unknown option $1"
            usage
            ;;
        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done

# Restore positional parameters
set -- "${POSITIONAL[@]}"

# Check if input file is provided
if [[ -z "$input_file" ]]; then
    echo "Error: No input file provided."
    usage
fi

# Verify if input file exists
if [[ ! -f "$input_file" ]]; then
    echo "Error: Input file not found: $input_file"
    exit 1
fi

# Set default output file if not provided
if [[ -z "$output_file" ]]; then
    base_name=$(basename "$input_file")
    base_name="${base_name%%.*}"
    output_file="../tmp/07.HGID/07.${base_name}.tsv"
fi

# Create output directory if it doesn't exist
mkdir -p "$(dirname "$output_file")"

# Define script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Define temporary directory
TMP_DIR="$SCRIPT_DIR/../tmp/07.HGID"
mkdir -p "$TMP_DIR"

# Temporary file names based on the input file name
base_name=$(basename "$input_file")
base_name="${base_name%%.*}"
tmp1="$TMP_DIR/tmp1_${base_name}.tsv"
tmp2="$TMP_DIR/tmp2_${base_name}.tsv"
tmp3="$TMP_DIR/tmp3_${base_name}.tsv"
tmp4="$TMP_DIR/tmp4_${base_name}.tsv"

# Echo progress
echo "Processing file: $input_file"
echo "Temporary files will be stored in: $TMP_DIR"

# Check if the input file is gzipped
if file "$input_file" | grep -q 'gzip compressed'; then
    echo "Unzipping and filtering VCF.gz file..."
    zcat "$input_file" | grep -v '^##' > "$tmp1"
else
    echo "Filtering header of VCF file..."
    grep -v '^##' "$input_file" > "$tmp1"
fi

# Check if required scripts are present
REQUIRED_SCRIPTS=("00.filter-vcf.py" "07.1.spliceai.py" "07.2.spliceai.py" "07.3.spliceai.py")
for script in "${REQUIRED_SCRIPTS[@]}"; do
    if [[ ! -e "$SCRIPT_DIR/$script" ]]; then
        echo "Error: Script not found: $script"
        exit 1
    fi
done

# Check if intervals file exists
if [[ ! -f "$intervals" ]]; then
    echo "Error: Intervals file not found: $intervals"
    exit 1
fi

#pipeline
echo "Filtering VCF data..."
python "$SCRIPT_DIR/00.filter-vcf.py" "$tmp1" "$intervals" > "$tmp2"
echo "Running 07.1.spliceai.py..."
python "$SCRIPT_DIR/07.1.spliceai.py" "$tmp2" "$tmp3"
echo "Running 07.2.spliceai.py..."
python "$SCRIPT_DIR/07.2.spliceai.py" "$tmp3" "$tmp4"
echo "Running 07.3.spliceai.py..."
python "$SCRIPT_DIR/07.3.spliceai.py" "$tmp4" "$output_file"

# Clean up temporary files
echo "Cleaning up temporary files..."
rm "$tmp1" "$tmp2" "$tmp3" "$tmp4"

# Final message
echo "Processing completed for $input_file."
echo "Output written to $output_file"
