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
    header_found = False
    try:
        with open_func(input_file_path, 'rt', encoding='utf-8-sig') as f_in:
            for line in f_in:
                clean_line = line.strip()
                if not clean_line:
                    continue
                # Skip VCF Metadata
                if clean_line.startswith('##'):
                    continue
                # Identify Header (Line starting with #CHROM or first non-comment line)
                if not header_found:
                    if clean_line.startswith('#'):
                        output_file.write(clean_line + '\n')
                        header_found = True
                        continue
                    else:
                        # For TSV files without # prefix
                        output_file.write(clean_line + '\n')
                        header_found = True
                        continue
                try:
                    parts = clean_line.split(sep)
                    if len(parts) < 2:
                        continue   
                    chrom, position = parts[0], parts[1]
                    # Attempt to parse position
                    pos_int = int(position) + 1
                    chrom_key = chrom[3:] if chrom.startswith('chr') else chrom
                    if chrom_key in bed_by_chrom:
                        regions = bed_by_chrom[chrom_key]
                        # Use bisect to find overlapping regions
                        index = bisect.bisect_right(regions, (pos_int, float('inf')))
                        if index > 0:
                            start, end = regions[index - 1]
                            if start - 4 <= pos_int <= end + 4:
                                output_file.write(clean_line + '\n')      
                except (ValueError, IndexError):
                    # Only log if it's not a header-looking line we missed
                    if not clean_line.startswith('#'):
                        logging.warning(f"Skipping malformed data line: {clean_line[:50]}...")
    except Exception as e:
        logging.error(f"Error reading input file: {e}")
        sys.exit(1)


def filter_input(input_file_path, bed_file_path, output_file_path=None):
    """
    Filters the input variant file based on bed regions, writing directly to output.

    Parameters:
    - input_file_path: Path to the variant file.
    - bed_file_path: Path to the bed data file.
    - output_file_path: Optional path to save the filtered output.
    """
    if input_file_path.endswith('.csv'):
        sep = ','
    else:
        sep = '\t'
    bed_by_chrom = read_bed_file(bed_file_path)
    try:
        if output_file_path:
            with open(output_file_path, 'w', encoding='utf-8') as f_out:
                batch_process_file(input_file_path, bed_by_chrom, sep, f_out)
        else:
            batch_process_file(input_file_path, bed_by_chrom, sep, sys.stdout)
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

