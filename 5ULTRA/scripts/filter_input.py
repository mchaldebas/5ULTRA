import bisect
import gzip

def read_bed_file(bed_file_path):
    """
    Reads the bed file and organizes regions by chromosome.
    """
    bed_by_chrom = {}
    with open(bed_file_path, 'r') as f:
        for line in f:
            chrom, start, end, *_ = line.strip().split('\t')
            chrom_key = chrom[3:] if chrom.startswith('chr') else chrom
            bed_by_chrom.setdefault(chrom_key, []).append((int(start), int(end)))
    for chrom in bed_by_chrom:
        bed_by_chrom[chrom].sort()
    return bed_by_chrom

def batch_process_file(input_file_path, bed_by_chrom, sep):
    """
    Processes the input variant file and filters variants based on bed regions.
    Returns a list of matched lines including the header.
    """
    matched_lines = []
    open_func = gzip.open if input_file_path.endswith('.gz') else open
    header = None
    with open_func(input_file_path, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if not header:
                header = line.strip()
                matched_lines.append(header)
                continue
            try:
                parts = line.rstrip('\n').split(sep, -1)
                chrom, position = parts[:2]
                position = int(position) + 1
                chrom_key = chrom[3:] if chrom.startswith('chr') else chrom
                if chrom_key in bed_by_chrom:
                    regions = bed_by_chrom[chrom_key]
                    index = bisect.bisect_right(regions, (position, float('inf')))
                    if index and regions[index - 1][0] <= position <= regions[index - 1][1]:
                        matched_lines.append(sep.join(parts))
            except:
                continue
    return matched_lines

def filter_input(input_file_path, bed_file_path, output_file_path=None):
    """
    Filters the input variant file based on bed regions.

    Parameters:
    - input_file_path: Path to the variant file.
    - bed_file_path: Path to the bed data file.
    - output_file_path: Optional path to save the filtered output.

    Returns:
    - List of matched lines if output_file_path is not provided.
    """
    sep = ',' if input_file_path.endswith('.csv') else '\t'
    bed_by_chrom = read_bed_file(bed_file_path)
    matched_lines = batch_process_file(input_file_path, bed_by_chrom, sep)
    if output_file_path:
        with open(output_file_path, 'w') as f_out:
            for line in matched_lines:
                f_out.write(line + '\n')
    else:
        return matched_lines

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Filter variants based on bed regions.')
    parser.add_argument('input_file_path', type=str, help='Path to the variant file.')
    parser.add_argument('bed_file_path', type=str, help='Path to the bed data file.')
    parser.add_argument('-o', '--output', type=str, help='Output file path (optional).')
    args = parser.parse_args()
    filter_input(args.input_file_path, args.bed_file_path, args.output)

if __name__ == "__main__":
    main()
