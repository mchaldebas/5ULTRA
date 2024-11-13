#!/bin/bash

set -e  # Exit immediately if a command exits with a non-zero status

# Function to print error messages and exit
error_exit() {
    echo "$1" >&2
    exit 1
}

# Function to calculate and print execution time
print_execution_time() {
    local description=$1
    local start_time=$2
    local end_time=$3
    echo -e "$description execution time:\t $(( end_time - start_time )) seconds"
}

# Trap to clean up temporary files on exit
cleanup() {
    rm -rf "$TMP_DIR"
}
trap cleanup EXIT

# Record start time
START_TIME=$(date +%s)

# Check if at least one argument is provided
if [ $# -lt 1 ]; then
    error_exit "Usage: $0 <input_file> [output_file]"
fi

# Set input and output files
input_file="$1"
if [ $# -ge 2 ]; then
    output_file="$2"
else
    output_file="${input_file%.*}.tsv"
fi
output="${output_file%.*}"

# Check if input file exists
if [ ! -f "$input_file" ]; then
    error_exit "Input file '$input_file' not found"
fi

# Create temporary directory
TMP_DIR=$(mktemp -d -t pipeline_tmp.XXXXXX)

# Unzip input file if necessary
step_start_time=$(date +%s)
if [[ "$input_file" == *.gz ]]; then
    unzipped_file="$TMP_DIR/$(basename "${input_file%.gz}")"
    gunzip -c "$input_file" > "$unzipped_file"
    input_file="$unzipped_file"
    step_end_time=$(date +%s)
    print_execution_time "Unzipping input file" $step_start_time $step_end_time
fi

# Check if required scripts and files exist
[ -f "./scripts/Filter-input.py" ] || error_exit "./scripts/Filter-input.py not found"
[ -f "./scripts/Detection.py" ] || error_exit "./scripts/Detection.py not found"
[ -f "./scripts/Score.py" ] || error_exit "./scripts/Score.py not found"
[ -f "./data/5UTRs.intervals.bed" ] || error_exit "Database file './Databases/5UTRs.intervals.bed' not found"

# Filter VCF/TSV data for 5'UTRs
step_start_time=$(date +%s)
filtered_output="$TMP_DIR/5UTR.${output##*/}.tsv"
python3 ./scripts/Filter-input.py "$input_file" "./data/5UTRs.intervals.bed" > "$filtered_output"
step_end_time=$(date +%s)
print_execution_time "5'UTR filtering" $step_start_time $step_end_time

# Run Detection
step_start_time=$(date +%s)
detection_output="$TMP_DIR/Detection.5UTR.${output##*/}.tsv"
python3 ./scripts/Detection.py "$filtered_output" "$detection_output"
step_end_time=$(date +%s)
print_execution_time "5'UTR detection" $step_start_time $step_end_time

# Run Scoring
step_start_time=$(date +%s)
python3 ./scripts/Score.py "$detection_output" "$output_file"
step_end_time=$(date +%s)
print_execution_time "Scoring" $step_start_time $step_end_time

# Calculate and print total execution time
END_TIME=$(date +%s)
execution_time=$(( END_TIME - START_TIME ))
echo -e "Total execution time:\t ${execution_time} seconds"
