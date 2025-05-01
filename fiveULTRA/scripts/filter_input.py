import bisect
import gzip
import logging
import argparse
import sys

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def read_bed_file(bed_file_path):
    """
    Reads the bed file and organizes regions by chromosome.
    """
    bed_by_chrom = {}
    try:
        with open(bed_file_path, 'r') as f:
            for line in f:
                if line.startswith("#"):
                    continue # Skip header lines in bed file
                chrom, start, end, *_ = line.strip().split('\t')
                chrom_key = chrom[3:] if chrom.startswith('chr') else chrom
                bed_by_chrom.setdefault(chrom_key, []).append((int(start), int(end)))
    except FileNotFoundError:
        logging.error(f"BED file not found: {bed_file_path}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error reading BED file: {e}")
        sys.exit(1)

    for chrom in bed_by_chrom:
        bed_by_chrom[chrom].sort()
    return bed_by_chrom


def batch_process_file(input_file_path, bed_by_chrom, sep, output_file):
    """
    Processes the input variant file and filters variants based on bed regions.
    Writes matched lines directly to the output file.
    """
    open_func = gzip.open if input_file_path.endswith('.gz') else open
    header = None
    try:
        with open_func(input_file_path, 'rt') as f_in:
            for line in f_in:
                if line.startswith('##'):
                    continue
                if not header:
                    header = line.strip()
                    output_file.write(header + '\n')
                    continue
                try:
                    parts = line.rstrip('\n').split(sep, -1)
                    chrom, position = parts[:2]
                    position = int(position) + 1
                    chrom_key = chrom[3:] if chrom.startswith('chr') else chrom
                    if chrom_key in bed_by_chrom:
                        regions = bed_by_chrom[chrom_key]
                        index = bisect.bisect_right(regions, (position, float('inf')))
                        if index and regions[index - 1][0] -4 <= position <= regions[index - 1][1] +4:
                            output_file.write(sep.join(parts) + '\n')
                except ValueError:
                    logging.warning(f"Skipping line due to ValueError: {line.strip()}")
                except IndexError:
                    logging.warning(f"Skipping line due to IndexError: {line.strip()}")

    except FileNotFoundError:
        logging.error(f"Input file not found: {input_file_path}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error reading input file: {e}")
        sys.exit(1)

    logging.info(f"Finished processing: {input_file_path}")


def filter_input(input_file_path, bed_file_path, output_file_path=None):
    """
    Filters the input variant file based on bed regions, writing directly to output.

    Parameters:
    - input_file_path: Path to the variant file.
    - bed_file_path: Path to the bed data file.
    - output_file_path: Optional path to save the filtered output.
    """

    sep = ',' if input_file_path.endswith('.csv') else '\t'
    bed_by_chrom = read_bed_file(bed_file_path)

    try:
        if output_file_path:
            with open(output_file_path, 'w') as f_out:
                batch_process_file(input_file_path, bed_by_chrom, sep, f_out)
            logging.info(f"Filtered data written to: {output_file_path}")
        else:
            # If no output path is given, write to standard output
            batch_process_file(input_file_path, bed_by_chrom, sep, sys.stdout)
            logging.info("Filtered data written to standard output.")

    except Exception as e:
        logging.error(f"Error during filtering: {e}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description='Filter variants based on bed regions.')
    parser.add_argument('input_file_path', type=str, help='Path to the variant file.')
    parser.add_argument('bed_file_path', type=str, help='Path to the bed data file.')
    parser.add_argument('-o', '--output', type=str, help='Output file path (optional, writes to stdout if not provided).')
    args = parser.parse_args()

    filter_input(args.input_file_path, args.bed_file_path, args.output)


if __name__ == "__main__":
    main()

