#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import tempfile
import gzip
import shutil
import time

def error_exit(message):
    """Print an error message to stderr and exit."""
    print(f"Error: {message}", file=sys.stderr)
    sys.exit(1)

def print_execution_time(description, start_time, end_time):
    """Print the execution time for a given step."""
    elapsed = int(end_time - start_time)
    print(f"{description} execution time:\t {elapsed} seconds")

def unzip_file(input_path, output_path):
    """Unzip a .gz file."""
    try:
        with gzip.open(input_path, 'rb') as f_in, open(output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    except OSError as e:
        error_exit(f"Failed to unzip '{input_path}': {e}")

def check_file_exists(path, description="File"):
    """Check if a file exists."""
    if not os.path.isfile(path):
        error_exit(f"{description} '{path}' not found.")

def run_subprocess(command, description):
    """Run a subprocess command and handle errors."""
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        error_exit(f"{description} failed with error: {e}")
    except FileNotFoundError:
        error_exit(f"Command not found: {' '.join(command)}")

def main():
    # Initialize argument parser
    parser = argparse.ArgumentParser(
        description="Process VCF/TSV files with optional SpliceAI processing.",
        usage="%(prog)s [--splice] <input_file> [output_file]"
    )
    parser.add_argument('--splice', action='store_true', help='Enable SpliceAI processing')
    parser.add_argument('input_file', help='Input VCF/TSV file')
    parser.add_argument('output_file', nargs='?', help='Output TSV file (default: input_file with .tsv extension)')

    args = parser.parse_args()

    splice = args.splice
    input_file = args.input_file
    output_file = args.output_file

    # Assign default output file if not provided
    if not output_file:
        base, _ = os.path.splitext(input_file)
        output_file = f"{base}.tsv"

    # Check if input file exists
    check_file_exists(input_file, "Input file")

    # Record start time
    start_time_total = time.time()

    # Create temporary directory
    with tempfile.TemporaryDirectory(prefix="pipeline_tmp.") as tmp_dir:
        # Unzip input file if necessary
        if input_file.endswith('.gz'):
            start_time = time.time()
            unzipped_filename = os.path.basename(input_file[:-3])
            unzipped_path = os.path.join(tmp_dir, unzipped_filename)
            print(f"Unzipping input file '{input_file}'...")
            unzip_file(input_file, unzipped_path)
            input_file = unzipped_path
            end_time = time.time()
            print_execution_time("Unzipping input file", start_time, end_time)

        # Check required scripts and files
        required_scripts = [
            "./scripts/Filter-input.py",
            "./scripts/Score.py"
        ]
        for script in required_scripts:
            check_file_exists(script, "Required script")

        required_data_file = "./data/5UTRs.intervals.bed"
        check_file_exists(required_data_file, "Database file")

        # Filter VCF/TSV data for 5'UTRs
        start_time = time.time()
        filtered_output = os.path.join(tmp_dir, f"5UTR.{os.path.basename(output_file)}.tsv")
        filter_command = [
            sys.executable, "./scripts/Filter-input.py",
            input_file,
            required_data_file
        ]
        print("Filtering input file for 5'UTRs...")
        with open(filtered_output, 'w') as f_out:
            try:
                subprocess.run(filter_command, check=True, stdout=f_out)
            except subprocess.CalledProcessError as e:
                error_exit(f"Filter-input.py failed with error: {e}")
            except FileNotFoundError:
                error_exit(f"Script not found: {' '.join(filter_command)}")
        end_time = time.time()
        print_execution_time("5'UTR filtering", start_time, end_time)

        # Conditional SpliceAI Processing
        if splice:
            spliceai_script = "./scripts/Spliceai-Main.py"
            check_file_exists(spliceai_script, "Script")
            spliceai_output = os.path.join(tmp_dir, f"Spliceai.{os.path.basename(output_file)}.tsv")
            spliceai_command = [
                spliceai_script,
                "--input", filtered_output,
                "--output", spliceai_output
            ]
            start_time = time.time()
            print("Running Spliceai-Main.py on detection output...")
            run_subprocess(spliceai_command, "Spliceai-Main.py")
            end_time = time.time()
            print_execution_time("Spliceai-Main.py", start_time, end_time)
            scoring_input = spliceai_output
        else:
            detection_script = "./scripts/Detection.py"
            check_file_exists(detection_script, "Script")
            detection_output = os.path.join(tmp_dir, f"Detection.5UTR.{os.path.basename(output_file)}.tsv")
            detection_command = [
                sys.executable, detection_script,
                filtered_output,
                detection_output
            ]
            start_time = time.time()
            print("Running Detection.py on filtered output...")
            run_subprocess(detection_command, "Detection.py")
            end_time = time.time()
            print_execution_time("5'UTR detection", start_time, end_time)
            scoring_input = detection_output

        # Run Scoring
        score_script = "./scripts/Score.py"
        score_command = [
            sys.executable, score_script,
            scoring_input,
            output_file
        ]
        start_time = time.time()
        print("Running Score.py...")
        run_subprocess(score_command, "Score.py")
        end_time = time.time()
        print_execution_time("Scoring", start_time, end_time)

    # Calculate and print total execution time
    end_time_total = time.time()
    total_time = int(end_time_total - start_time_total)
    print(f"Total execution time:\t {total_time} seconds")

if __name__ == "__main__":
    main()
