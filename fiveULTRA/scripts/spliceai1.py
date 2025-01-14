import pysam
import os
import sys

def parse_spliceai_info(info_field):
    """
    Extracts the SpliceAI annotation from the INFO field of a VCF record.
    """
    for entry in info_field.split(';'):
        if entry.startswith('SpliceAI='):
            return entry[len('SpliceAI='):]
    return None

def process_variant_spliceai_1(chrom, pos, ref, alt, snv_vcf, indel_vcf):
    """
    Processes a single variant by querying SpliceAI annotations from VCF files.

    Parameters:
    - chrom: Chromosome (string).
    - pos: Position (int).
    - ref: Reference allele (string).
    - alt: Alternate allele (string).
    - snv_vcf: pysam.TabixFile object for SNV VCF.
    - indel_vcf: pysam.TabixFile object for indel VCF.

    Returns:
    - List of SpliceAI annotations.
    """
    spliceai_results = []

    def query_vcf(vcf_file):
        try:
            records = vcf_file.fetch(chrom, pos - 1, pos)
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
        except (ValueError, IndexError, pysam.TabixError):
            pass  # No records found or invalid data

    # Query SNV VCF
    if len(ref) == 1 and len(alt) == 1:
        query_vcf(snv_vcf)
    # If not found, query indel VCF
    else:
        query_vcf(indel_vcf)

    return spliceai_results

def process_spliceai_1(input_file, output_file, data_dir='~/.5ULTRA/data', cutoff=0.2):
    """
    Processes an input file to add SpliceAI annotations.

    Parameters:
    - input_file: Path to the input file.
    - output_file: Path to the output file.
    - data_dir: Directory containing SpliceAI data files.
    - cutoff: Score cutoff for SpliceAI annotations.
    """
    # Verify if input file exists
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Input file '{input_file}' not found.")

    # Paths to SpliceAI VCF files
    snv_vcf_path = os.path.join(os.path.expanduser(data_dir), "spliceai50.5UTRs.raw.snvs.hg38.vcf.gz")
    indel_vcf_path = os.path.join(os.path.expanduser(data_dir), "spliceai50.5UTRs.raw.indels.hg38.vcf.gz")

    if not os.path.isfile(snv_vcf_path):
        raise FileNotFoundError(f"SNV VCF file '{snv_vcf_path}' not found.")
    if not os.path.isfile(indel_vcf_path):
        raise FileNotFoundError(f"Indel VCF file '{indel_vcf_path}' not found.")

    # Open VCF files using pysam
    snv_vcf = pysam.TabixFile(snv_vcf_path)
    indel_vcf = pysam.TabixFile(indel_vcf_path)

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Write the header line to the output file
        header = infile.readline().strip('\r\n')
        outfile.write(f"{header}\tSpliceAI\n")

        # Process each line in the input file
        for line in infile:
            line = line.strip('\r\n')
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) < 5:
                continue  # Skip incomplete lines

            chrom = fields[0].lstrip('chr')
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]

            # Process variant
            spliceai_annotations = process_variant_spliceai_1(chrom, pos, ref, alt, snv_vcf, indel_vcf)

            # Process each SpliceAI annotation
            for annotation in spliceai_annotations:
                entries = annotation.split(',')
                for entry in entries:
                    try:
                        allele, gene, score1, score2, score3, score4, pos1, pos2, pos3, pos4 = entry.split('|')
                        scores = [float(score1), float(score2), float(score3), float(score4)]
                        positions = [pos1, pos2, pos3, pos4]
                        if any(score > cutoff for score in scores):
                            formatted_scores = '|'.join([gene] + [f"{score:.2f}" for score in scores] + positions)
                            outfile.write(f"{line}\t{formatted_scores}\n")
                    except ValueError:
                        continue  # Skip entries with invalid format

    # Close VCF files
    snv_vcf.close()
    indel_vcf.close()

# Optional main function to allow script execution directly
def main():
    import argparse
    parser = argparse.ArgumentParser(description='Process SpliceAI annotations.')
    parser.add_argument('input_file', type=str, help='Path to the input file.')
    parser.add_argument('output_file', type=str, help='Path to the output file.')
    parser.add_argument('--data-dir', type=str, default='~/.5ULTRA/data', help='Path to the data directory.')
    parser.add_argument('--cutoff', type=float, default=0.2, help='Score cutoff for SpliceAI annotations.')
    args = parser.parse_args()
    try:
        process_spliceai_1(args.input_file, args.output_file, data_dir=args.data_dir, cutoff=args.cutoff)
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
