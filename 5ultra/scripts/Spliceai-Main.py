#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import tempfile
import time

def error_exit(message):
    """Print an error message to stderr and exit."""
    print(f"Error: {message}", file=sys.stderr)
    sys.exit(1)

def print_execution_time(description, start_time, end_time):
    """Print the execution time for a given step."""
    elapsed = int(end_time - start_time)
    print(f"{description} execution time:\t {elapsed} seconds")

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
        description="Process VCF/TSV files with SpliceAI processing.",
        usage="%(prog)s <input_file> <output_file>"
    )
    parser.add_argument('input_file', help='Input VCF/TSV file')
    parser.add_argument('output_file', help='Output TSV file')
    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file
    # Check if input file exists
    check_file_exists(input_file, "Input file")
    # Record start time
    start_time_total = time.time()
    # Create temporary directory
    with tempfile.TemporaryDirectory(prefix="pipeline_tmp.") as tmp_dir:
        # Check required scripts and files
        required_scripts = [
            "./scripts/Spliceai-1.py",
            "./scripts/Spliceai-2.py",
            "./scripts/Spliceai-3.py"
        ]
        for script in required_scripts:
            check_file_exists(script, "Required script")

        spliceai_script = "./scripts/Spliceai-1.py"
        check_file_exists(spliceai_script, "Script")
        spliceai_output = os.path.join(tmp_dir, f"Spliceai-1.{os.path.basename(output_file)}")
        spliceai_command = [ sys.executable,
            spliceai_script,
            input_file,
            spliceai_output
        ]
        start_time = time.time()
        print("Running Spliceai-1.py on filtering output...")
        run_subprocess(spliceai_command, "Spliceai-1.py")
        end_time = time.time()
        print_execution_time("Spliceai-1.py", start_time, end_time)
        next_input = spliceai_output

        spliceai_script = "./scripts/Spliceai-2.py"
        check_file_exists(spliceai_script, "Script")
        spliceai_output = os.path.join(tmp_dir, f"Spliceai-2.{os.path.basename(output_file)}")
        spliceai_command = [ sys.executable,
            spliceai_script,
            next_input,
            spliceai_output
        ]
        start_time = time.time()
        print("Running Spliceai-2.py on Spliceai-1 output...")
        run_subprocess(spliceai_command, "Spliceai-2.py")
        end_time = time.time()
        print_execution_time("Spliceai-2.py", start_time, end_time)
        next_input = spliceai_output

        spliceai_script = "./scripts/Spliceai-3.py"
        check_file_exists(spliceai_script, "Script")
        spliceai_output = output_file
        spliceai_command = [ sys.executable,
            spliceai_script,
            next_input,
            spliceai_output
        ]
        start_time = time.time()
        print("Running Spliceai-3.py on Spliceai-2 output...")
        run_subprocess(spliceai_command, "Spliceai-3.py")
        end_time = time.time()
        print_execution_time("Spliceai-3.py", start_time, end_time)

    # Calculate and print total execution time
    end_time_total = time.time()
    total_time = int(end_time_total - start_time_total)
    print(f"Total execution time:\t {total_time} seconds")

if __name__ == "__main__":
    main()
