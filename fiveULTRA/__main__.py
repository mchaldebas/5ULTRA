import sys
import argparse
import logging
import os
import tempfile
import gzip
import shutil
import time

from .scripts.filter_input import filter_input
from .scripts.detection import process_variants, load_vcf_data, load_tsv_data
from .scripts.score import score_variants
from .scripts.spliceai1 import process_spliceai_1
from .scripts.spliceai2 import process_variants_spliceai_2
from .scripts.spliceai3 import process_variants_spliceai_3

import argparse

def get_options():
    parser = argparse.ArgumentParser(
        prog="5ULTRA",
        description=(
            "5ULTRA: Computational pipeline designed to annotate and score genetic variants located in the 5′UTR"
        ),
        epilog="Example: 5ULTRA -I input.vcf -O output.tsv --splice --full",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # Input/Output arguments
    io_group = parser.add_argument_group("Input/Output options")
    io_group.add_argument('-I', metavar='input', required=True,
                          help='Path to the input VCF/TSV file (required)')
    io_group.add_argument('-O', metavar='output',
                          help='Path to the output TSV file (optional)')
    # Processing options
    proc_group = parser.add_argument_group("Processing options")
    proc_group.add_argument('--splice', action='store_true',
                            help='Enable SpliceAI processing')
    proc_group.add_argument('--full', action='store_true',
                            help='Enable full annotation')
    proc_group.add_argument('--mane', action='store_true',
                            help='Focus on MANE (Matched Annotation from the NCBI and EMBL-EBI) transcripts')
    # Advanced options
    advanced_group = parser.add_argument_group("Advanced options")
    advanced_group.add_argument('--data-dir', metavar='data_path', type=str,
                                default='~/.5ULTRA/data',
                                help='Path to the data directory for SpliceAI resources')

    args = parser.parse_args()
    return args

def main():
    args = get_options()
    splice = args.splice
    input_file = args.I
    output_file = args.O
    data_dir = args.data_dir
    full_anno = args.full
    mane = args.mane
    cutoff = 0.2
    # Assign default output file if not provided
    if not output_file:
        base, _ = os.path.splitext(input_file)
        output_file = f"{base}.5ULTRA.tsv"

    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(message)s')

    # Check if input file exists
    if not os.path.isfile(input_file):
        logging.error(f"Input file '{input_file}' not found.")
        sys.exit(1)

    # Record start time
    start_time_total = time.time()

    # Create temporary directory
    with tempfile.TemporaryDirectory(prefix="pipeline_tmp.") as tmp_dir:
        # Unzip input file if necessary
        if input_file.endswith('.gz'):
            start_time = time.time()
            unzipped_filename = os.path.basename(input_file[:-3])
            unzipped_path = os.path.join(tmp_dir, unzipped_filename)
            logging.info(f"Unzipping input file '{input_file}'...")
            try:
                with gzip.open(input_file, 'rb') as f_in, open(unzipped_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                input_file = unzipped_path
            except OSError as e:
                logging.error(f"Failed to unzip '{input_file}': {e}")
                sys.exit(1)
            end_time = time.time()
            logging.info(f"Unzipping input file execution time:\t {int(end_time - start_time)} seconds")

        # Paths to scripts and data
        required_data_file = os.path.join(os.path.expanduser(data_dir), '5UTRs.intervals.bed')
        if not os.path.isfile(required_data_file):
            logging.error(f"Database file '{required_data_file}' not found.")
            sys.exit(1)

        # Filter VCF/TSV data for 5'UTRs
        start_time = time.time()
        filtered_output = os.path.join(tmp_dir, f"5UTR.{os.path.basename(output_file)}.tsv")
        logging.info("Filtering input file for 5'UTRs...")
        try:
            filter_input(input_file, required_data_file, filtered_output)
        except Exception as e:
            logging.error(f"Filter-input failed with error: {e}")
            sys.exit(1)

        end_time = time.time()
        logging.info(f"5'UTR filtering execution time:\t {int(end_time - start_time)} seconds")

        # Conditional SpliceAI Processing
        if splice:
            # spliceai Detection processing
            output_file = f"{base}.5ULTRA_splice.tsv"
            splice_1_output = os.path.join(tmp_dir, f"splice1.5UTR.{os.path.basename(output_file)}")
            splice_2_output = os.path.join(tmp_dir, f"splice2.5UTR.{os.path.basename(output_file)}")
            splice_3_output = os.path.join(tmp_dir, f"splice3.5UTR.{os.path.basename(output_file)}")
            start_time = time.time()
            logging.info("Running splice detection on filtered output...")
            try:
                process_spliceai_1(filtered_output, splice_1_output, data_dir, cutoff)
                process_variants_spliceai_2(splice_1_output, splice_2_output, data_dir, cutoff)
                process_variants_spliceai_3(splice_2_output, splice_3_output, data_dir)
            except Exception as e:
                logging.error(f"splice Detection failed with error: {e}")
                sys.exit(1)
            end_time = time.time()
            logging.info(f"5'UTR splice detection execution time:\t {int(end_time - start_time)} seconds")
            scoring_input = splice_3_output
        else:
            # Detection processing
            detection_output = os.path.join(tmp_dir, f"Detection.5UTR.{os.path.basename(output_file)}")
            start_time = time.time()
            logging.info("Running detection on filtered output...")
            try:
                if filtered_output.endswith('.vcf'):
                    variants = load_vcf_data(filtered_output)
                else:
                    variants = load_tsv_data(filtered_output)
                process_variants(variants, detection_output, data_dir)
            except Exception as e:
                logging.error(f"Detection failed with error: {e}")
                sys.exit(1)
            end_time = time.time()
            logging.info(f"5'UTR detection execution time:\t {int(end_time - start_time)} seconds")
            scoring_input = detection_output

        # Run Scoring
        start_time = time.time()
        logging.info("Running scoring...")
        try:
            res = score_variants(scoring_input, output_file, data_dir, full_anno, mane)
            if not res:
                sys.exit(1)
        except Exception as e:
            logging.error(f"Scoring failed with error: {e}")
            sys.exit(1)

        end_time = time.time()
        logging.info(f"Scoring execution time:\t {int(end_time - start_time)} seconds")

    # Calculate and print total execution time
    end_time_total = time.time()
    total_time = int(end_time_total - start_time_total)
    logging.info(f"Total execution time:\t {total_time} seconds")
    logging.info(f"Results available:\t {output_file}")

if __name__ == '__main__':
    main()
