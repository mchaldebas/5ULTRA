#!/bin/bash

set -e  # Exit immediately if a command exits with a non-zero status

# Function to print error messages and exit
error_exit() {
    echo "$1" >&2
    exit 1
}

# Function to calculate and print execution time
print_execution_time() {
    local description="$1"
    local start_time="$2"
    local end_time="$3"
    echo -e "$description execution time:\t $(( end_time - start_time )) seconds"
}

# Function to display usage instructions
usage() {
    echo "Usage: $0 [--splice] <input_file> [output_file]"
    echo
    echo "Options:"
    echo "  --splice          Enable SpliceAI processing"
    echo "  -h, --help        Display this help message"
    exit 1
}

# Initialize variables
splice=false
input_file=""
output_file=""

# Parse command-line arguments
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --splice)
            splice=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        --)
            shift
            while [[ $# -gt 0 ]]; do
                POSITIONAL+=("$1")
                shift
            done
            ;;
        -*)
            echo "Unknown option: $1"
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

# Assign input and output files
if [[ $# -lt 1 ]]; then
    echo "Error: Missing input file."
    usage
fi

input_file="$1"
if [[ $# -ge 2 ]]; then
    output_file="$2"
else
    output_file="${input_file%.*}.tsv"
fi

# Check if input file exists
if [[ ! -f "$input_file" ]]; then
    error_exit "Input file '$input_file' not found."
fi

# Record start time
START_TIME=$(date +%s)

# Create temporary directory
TMP_DIR=$(mktemp -d -t pipeline_tmp.XXXXXX)

# Ensure temporary directory is cleaned up on exit
cleanup() {
    rm -rf "$TMP_DIR"
}
trap cleanup EXIT

# Unzip input file if necessary
step_start_time=$(date +%s)
if [[ "$input_file" == *.gz ]]; then
    unzipped_file="$TMP_DIR/$(basename "${input_file%.gz}")"
    gunzip -c "$input_file" > "$unzipped_file"
    input_file="$unzipped_file"
    step_end_time=$(date +%s)
    print_execution_time "Unzipping input file" "$step_start_time" "$step_end_time"
fi

# Check if required scripts and files exist
REQUIRED_SCRIPTS=(
    "./scripts/Filter-input.py"
    "./scripts/Score.py"
)
for script in "${REQUIRED_SCRIPTS[@]}"; do
    [[ -f "$script" ]] || error_exit "Required script not found: $script"
done

[[ -f "./data/5UTRs.intervals.bed" ]] || error_exit "Database file './data/5UTRs.intervals.bed' not found."

# Filter VCF/TSV data for 5'UTRs
step_start_time=$(date +%s)
filtered_output="$TMP_DIR/5UTR.${output_file##*/}.tsv"
python3 ./scripts/Filter-input.py "$input_file" "./data/5UTRs.intervals.bed" > "$filtered_output"
step_end_time=$(date +%s)
print_execution_time "5'UTR filtering" "$step_start_time" "$step_end_time"

# Conditional SpliceAI Processing
if $splice; then
    # Define the path to Spliceai-Main.sh
    spliceai_script="./scripts/Spliceai-Main.sh"
    # Check if Spliceai-Main.sh exists
    [[ -f "$spliceai_script" ]] || error_exit "Script not found: $spliceai_script"
    # Define the output file for Spliceai-Main.sh
    spliceai_output="$TMP_DIR/Spliceai.${output_file##*/}.tsv"
    # Run Spliceai-Main.sh
    step_start_time=$(date +%s)
    echo "Running Spliceai-Main.sh on detection output..."
    bash "$spliceai_script" --input "$filtered_output" --output "$spliceai_output"
    step_end_time=$(date +%s)
    print_execution_time "Spliceai-Main.sh" "$step_start_time" "$step_end_time"
    # Set the input for Scoring
    scoring_input="$spliceai_output"
else
    Detection_script="./scripts/Detection.py"
    # Check if Spliceai-Main.sh exists
    [[ -f "$Detection_script" ]] || error_exit "Script not found: $Detection_script"
    # Run Detection
    step_start_time=$(date +%s)
    detection_output="$TMP_DIR/Detection.5UTR.${output_file##*/}.tsv"
    python3 "$Detection_script" "$filtered_output" "$detection_output"
    step_end_time=$(date +%s)
    print_execution_time "5'UTR detection" "$step_start_time" "$step_end_time"
    # Set the input for Scoring
    scoring_input="$detection_output"
fi

# Run Scoring
step_start_time=$(date +%s)
python3 ./scripts/Score.py "$scoring_input" "$output_file"
step_end_time=$(date +%s)
print_execution_time "Scoring" "$step_start_time" "$step_end_time"

# Calculate and print total execution time
END_TIME=$(date +%s)
execution_time=$(( END_TIME - START_TIME ))
echo -e "Total execution time:\t ${execution_time} seconds"
