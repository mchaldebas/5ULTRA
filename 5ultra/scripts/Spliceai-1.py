#!/usr/bin/env python3

import sys
import pysam

def usage():
    print(f"Usage: {sys.argv[0]} input_file output_file")
    sys.exit(1)

if len(sys.argv) != 3:
    usage()

input_file = sys.argv[1]
output_file = sys.argv[2]
cutoff = 0.2

# Verify if input file exists
try:
    infile = open(input_file, 'r')
except FileNotFoundError:
    print("Input file not found!")
    sys.exit(1)

# Open output file
outfile = open(output_file, 'w')

# Write the header line to the output file
header = infile.readline().strip('\r\n')
outfile.write(f"{header}\tSpliceAI\n")

# Paths to SpliceAI VCF files (update these paths accordingly)
snv_vcf_path = "./data/spliceai_scores.raw.snv.hg38.vcf.gz"
indel_vcf_path = "./data/spliceai_scores.raw.indel.hg38.vcf.gz"

# Open VCF files using pysam
snv_vcf = pysam.TabixFile(snv_vcf_path)
indel_vcf = pysam.TabixFile(indel_vcf_path)

def parse_spliceai_info(info_field):
    # Extract the SpliceAI annotation
    for entry in info_field.split(';'):
        if entry.startswith('SpliceAI='):
            return entry[len('SpliceAI='):]
    return None

def process_variant(chr, pos, ref, alt):
    spliceai_results = []

    # Function to query VCF and extract SpliceAI annotations
    def query_vcf(vcf_file):
        try:
            records = vcf_file.fetch(chr, pos-1, pos)
            for record in records:
                fields = record.strip().split('\t')
                vcf_pos = int(fields[1])
                vcf_ref = fields[3]
                vcf_alt = fields[4]
                info = fields[7]
                if vcf_pos == pos and vcf_ref == ref and vcf_alt == alt:
                    spliceai_annotation = parse_spliceai_info(info)
                    if spliceai_annotation:
                        spliceai_results.append(spliceai_annotation)
        except ValueError:
            pass  # No records found

    # Query SNV VCF
    query_vcf(snv_vcf)
    # If not found, query indel VCF
    if not spliceai_results:
        query_vcf(indel_vcf)

    return spliceai_results

# Process each line in the input file
for line in infile:
    line = line.strip('\r\n')
    if not line:
        continue
    fields = line.split('\t')
    if len(fields) < 5:
        continue  # Skip incomplete lines

    chr = fields[0].lstrip('chr')
    pos = int(fields[1])
    ref = fields[3]
    alt = fields[4]

    # Process variant
    spliceai_annotations = process_variant(chr, pos, ref, alt)

    # Process each SpliceAI annotation
    for annotation in spliceai_annotations:
        entries = annotation.split(',')
        for entry in entries:
            allele, gene, score1, score2, score3, score4, pos1, pos2, pos3, pos4 = entry.split('|')
            scores = [float(score1), float(score2), float(score3), float(score4)]
            positions = [pos1, pos2, pos3, pos4]
            if any(score > cutoff for score in scores):
                formatted_scores = '|'.join([gene] + [f"{score:.2f}" for score in scores] + positions)
                outfile.write(f"{line}\t{formatted_scores}\n")

# Close files
infile.close()
outfile.close()
