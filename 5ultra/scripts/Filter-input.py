#!/usr/bin/env python3

import argparse
import bisect
import gzip

def read_introns_file(introns_file_path):
    introns_by_chrom = {}
    with open(introns_file_path, 'r') as f:
        for line in f:
            chrom, start, end, *_ = line.strip().split('\t')
            chrom_key = chrom[3:] if chrom.startswith('chr') else chrom
            introns_by_chrom.setdefault(chrom_key, []).append((int(start), int(end)))
    for chrom in introns_by_chrom:
        introns_by_chrom[chrom].sort()
    return introns_by_chrom

def batch_process_file(input_file_path, introns_by_chrom, sep):
    matched_lines = []
    open_func = gzip.open if input_file_path.endswith('.gz') else open
    header = None
    with open_func(input_file_path, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if not header:
                header = line
                print('\t'.join(line.strip().split(sep)))
                continue
            try:
                parts = line.rstrip('\n').split(sep, -1)
                chrom, position = parts[:2]
                position = int(position) +1
                chrom_key = chrom[3:] if chrom.startswith('chr') else chrom
                if chrom_key in introns_by_chrom:
                    introns = introns_by_chrom[chrom_key]
                    index = bisect.bisect_right(introns, (position, float('inf')))
                    if index and introns[index - 1][0] <= position <= introns[index - 1][1]:
                        print(sep.join(parts))
            except: continue

def main():
    parser = argparse.ArgumentParser(description='Filter variants based on intron regions.')
    parser.add_argument('input_file_path', type=str, help='Path to the variant file.')
    parser.add_argument('introns_file_path', type=str, help='Path to the intron data file.')
    args = parser.parse_args()
    introns_by_chrom = read_introns_file(args.introns_file_path)
    sep = ',' if args.input_file_path.endswith('.csv') else '\t'
    batch_process_file(args.input_file_path, introns_by_chrom, sep)

if __name__ == "__main__":
    main()

