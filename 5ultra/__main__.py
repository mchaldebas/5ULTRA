import sys
import argparse
import logging
import os
import subprocess
import tempfile
import gzip
import shutil
import time

def get_options():
    parser = argparse.ArgumentParser(
        description="Process VCF/TSV files with optional SpliceAI processing."
    )
    parser.add_argument('-I', metavar='input', nargs='?', required=True,
                        help='Path to the input VCF/TSV file')
    parser.add_argument('-O', metavar='output', nargs='?', help='Path to the output TSV file')
    parser.add_argument('--splice', action='store_true', help='Enable SpliceAI processing')
    args = parser.parse_args()
    return args

def main():
    args = get_options()
    splice = args.splice
    input_file = args.I
    output_file = args.O

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

        # Check required scripts and files
        required_scripts = [
            "./scripts/Filter-input.py",
            "./scripts/Score.py"
        ]
        for script in required_scripts:
            if not os.path.isfile(script):
                logging.error(f"Required script '{script}' not found.")
                sys.exit(1)

        required_data_file = "./data/5UTRs.intervals.bed"
        if not os.path.isfile(required_data_file):
            logging.error(f"Database file '{required_data_file}' not found.")
            sys.exit(1)

        # Filter VCF/TSV data for 5'UTRs
        start_time = time.time()
        filtered_output = os.path.join(tmp_dir, f"5UTR.{os.path.basename(output_file)}.tsv")
        filter_command = [
            sys.executable, "./scripts/Filter-input.py",
            input_file,
            required_data_file
        ]
        logging.info("Filtering input file for 5'UTRs...")
        with open(filtered_output, 'w') as f_out:
            try:
                subprocess.run(filter_command, check=True, stdout=f_out)
            except subprocess.CalledProcessError as e:
                logging.error(f"Filter-input.py failed with error: {e}")
                sys.exit(1)
            except FileNotFoundError:
                logging.error(f"Script not found: {' '.join(filter_command)}")
                sys.exit(1)
        end_time = time.time()
        logging.info(f"5'UTR filtering execution time:\t {int(end_time - start_time)} seconds")

        # Conditional SpliceAI Processing
        if splice:
            spliceai_script = "./scripts/Spliceai-Main.py"
            if not os.path.isfile(spliceai_script):
                logging.error(f"Script '{spliceai_script}' not found.")
                sys.exit(1)
            spliceai_output = os.path.join(tmp_dir, f"Spliceai.{os.path.basename(output_file)}")
            spliceai_command = [sys.executable,
                spliceai_script,
                filtered_output,
                spliceai_output
            ]
            start_time = time.time()
            logging.info("Running Spliceai-Main.py on detection output...")
            try:
                subprocess.run(spliceai_command, check=True)
            except subprocess.CalledProcessError as e:
                logging.error(f"Spliceai-Main.py failed with error: {e}")
                sys.exit(1)
            except FileNotFoundError:
                logging.error(f"Script not found: {' '.join(spliceai_command)}")
                sys.exit(1)
            end_time = time.time()
            logging.info(f"Spliceai-Main.py execution time:\t {int(end_time - start_time)} seconds")
            scoring_input = spliceai_output
        else:
            detection_script = "./scripts/Detection.py"
            if not os.path.isfile(detection_script):
                logging.error(f"Script '{detection_script}' not found.")
                sys.exit(1)
            detection_output = os.path.join(tmp_dir, f"Detection.5UTR.{os.path.basename(output_file)}")
            detection_command = [
                sys.executable, detection_script,
                filtered_output,
                detection_output
            ]
            start_time = time.time()
            logging.info("Running Detection.py on filtered output...")
            try:
                subprocess.run(detection_command, check=True)
            except subprocess.CalledProcessError as e:
                logging.error(f"Detection.py failed with error: {e}")
                sys.exit(1)
            except FileNotFoundError:
                logging.error(f"Script not found: {' '.join(detection_command)}")
                sys.exit(1)
            end_time = time.time()
            logging.info(f"5'UTR detection execution time:\t {int(end_time - start_time)} seconds")
            scoring_input = detection_output

        # Run Scoring
        score_script = "./scripts/Score.py"
        if not os.path.isfile(score_script):
            logging.error(f"Script '{score_script}' not found.")
            sys.exit(1)
        score_command = [
            sys.executable, score_script,
            scoring_input,
            output_file
        ]
        start_time = time.time()
        logging.info("Running Score.py...")
        try:
            subprocess.run(score_command, check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Score.py failed with error: {e}")
            sys.exit(1)
        except FileNotFoundError:
            logging.error(f"Script not found: {' '.join(score_command)}")
            sys.exit(1)
        end_time = time.time()
        logging.info(f"Scoring execution time:\t {int(end_time - start_time)} seconds")

    # Calculate and print total execution time
    end_time_total = time.time()
    total_time = int(end_time_total - start_time_total)
    logging.info(f"Total execution time:\t {total_time} seconds")

if __name__ == '__main__':
    main()
