#!/usr/bin/env python3

import argparse
import sys
import gzip
from pathlib import Path
import subprocess

def print_error(message):
    print(f"Error: {message}", file=sys.stderr)

def usage():
    print("Usage: script.py --input input_file [--output output_file]", file=sys.stderr)
    sys.exit(1)

def is_gzipped(file_path):
    with open(file_path, 'rb') as f:
        magic_number = f.read(2)
    return magic_number == b'\x1f\x8b'

def filter_vcf(input_file, tmp1, is_gz):
    print("Filtering VCF header...")
    try:
        if is_gz:
            with gzip.open(input_file, 'rt') as fin, open(tmp1, 'w') as fout:
                for line in fin:
                    if not line.startswith('##'):
                        fout.write(line)
        else:
            with open(input_file, 'r') as fin, open(tmp1, 'w') as fout:
                for line in fin:
                    if not line.startswith('##'):
                        fout.write(line)
    except Exception as e:
        print_error(f"Failed to filter VCF file: {e}")
        sys.exit(1)

def check_scripts(script_dir, required_scripts):
    missing_scripts = []
    for script in required_scripts:
        script_path = script_dir / script
        if not script_path.is_file():
            missing_scripts.append(script)
    if missing_scripts:
        print_error(f"Script(s) not found: {', '.join(missing_scripts)}")
        sys.exit(1)

def check_intervals(intervals_path):
    if not intervals_path.is_file():
        print_error(f"Intervals file not found: {intervals_path}")
        sys.exit(1)

def run_pipeline(script_dir, tmp_dir, tmp_files, output_file):
    try:
        # Define the sequence of scripts and their arguments
        pipeline_steps = [
            {
                'script': '00.filter-vcf.py',
                'args': [str(tmp_files['tmp1']), str(tmp_files['intervals'])],
                'output': tmp_files['tmp2']
            },
            {
                'script': '07.1.spliceai.py',
                'args': [str(tmp_files['tmp2']), str(tmp_files['tmp3'])],
                'output': tmp_files['tmp3']
            },
            {
                'script': '07.2.spliceai.py',
                'args': [str(tmp_files['tmp3']), str(tmp_files['tmp4'])],
                'output': tmp_files['tmp4']
            },
            {
                'script': '07.3.spliceai.py',
                'args': [str(tmp_files['tmp4']), str(output_file)],
                'output': output_file
            }
        ]

        for step in pipeline_steps:
            script_path = script_dir / step['script']
            print(f"Running {step['script']}...")
            subprocess.run(
                ['python', str(script_path)] + step['args'],
                check=True
            )
    except subprocess.CalledProcessError as e:
        print_error(f"Pipeline step failed: {e}")
        sys.exit(1)
    except Exception as e:
        print_error(f"Unexpected error during pipeline: {e}")
        sys.exit(1)

def main():
    # Default values
    default_intervals = Path("./data/5UTRs.intervals.bed")

    # Argument parsing
    parser = argparse.ArgumentParser(description="Process VCF files.")
    parser.add_argument('--input', required=True, help='Input VCF file')
    parser.add_argument('--output', help='Output TSV file')
    args = parser.parse_args()

    input_file = Path(args.input)
    output_file = Path(args.output) if args.output else None
    intervals = default_intervals

    # Check if input file is provided
    if not input_file:
        print_error("No input file provided.")
        usage()

    # Verify if input file exists
    if not input_file.is_file():
        print_error(f"Input file not found: {input_file}")
        sys.exit(1)

    # Set default output file if not provided
    if not output_file:
        base_name = input_file.stem
        output_file = Path("../tmp/07.HGID") / f"07.{base_name}.tsv"

    # Create output directory if it doesn't exist
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Define script directory
    try:
        script_dir = Path(__file__).resolve().parent
    except NameError:
        # If __file__ is not defined, use current working directory
        script_dir = Path.cwd()

    # Define temporary directory
    tmp_dir = script_dir.parent / "tmp" / "07.HGID"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # Temporary file names based on the input file name
    base_name = input_file.stem
    tmp1 = tmp_dir / f"tmp1_{base_name}.tsv"
    tmp2 = tmp_dir / f"tmp2_{base_name}.tsv"
    tmp3 = tmp_dir / f"tmp3_{base_name}.tsv"
    tmp4 = tmp_dir / f"tmp4_{base_name}.tsv"

    tmp_files = {
        'tmp1': tmp1,
        'tmp2': tmp2,
        'tmp3': tmp3,
        'tmp4': tmp4,
        'intervals': intervals
    }

    # Echo progress
    print(f"Processing file: {input_file}")
    print(f"Temporary files will be stored in: {tmp_dir}")

    # Check if the input file is gzipped
    is_gz = is_gzipped(input_file)
    if is_gz:
        print("Unzipping and filtering VCF.gz file...")
    else:
        print("Filtering header of VCF file...")

    # Filter VCF and create tmp1
    filter_vcf(input_file, tmp1, is_gz)

    # Check if required scripts are present
    required_scripts = ["00.filter-vcf.py", "07.1.spliceai.py", "07.2.spliceai.py", "07.3.spliceai.py"]
    check_scripts(script_dir, required_scripts)

    # Check if intervals file exists
    check_intervals(intervals)

    # Run the pipeline
    run_pipeline(script_dir, tmp_dir, tmp_files, output_file)

    # Clean up temporary files
    print("Cleaning up temporary files...")
    try:
        tmp1.unlink()
        tmp2.unlink()
        tmp3.unlink()
        tmp4.unlink()
    except Exception as e:
        print_error(f"Failed to remove temporary files: {e}")
        # Not exiting since processing is complete

    # Final message
    print(f"Processing completed for {input_file}.")
    print(f"Output written to {output_file}")

if __name__ == "__main__":
    main()
