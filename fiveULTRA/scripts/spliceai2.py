import csv
import sys
import os
from collections import defaultdict
import ast

# Constants
BASE_PAIRING = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N', '*': '*'}

def rev_seq(fwd_seq):
    """Reverse complements a DNA sequence based on the BASE_PAIRING dictionary."""
    return ''.join(BASE_PAIRING.get(nuc, 'N') for nuc in reversed(fwd_seq))

def load_tsv_data(file_path):
    """Loads TSV data from a file and returns a list of rows."""
    data = []
    try:
        max_int = sys.maxsize
        while True:
            try:
                csv.field_size_limit(max_int)
                break
            except OverflowError:
                max_int = int(max_int / 10)
        with open(file_path, 'r') as file:
            csv_reader = csv.reader(file, delimiter='\t')
            data = list(csv_reader)
    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
    return data

def calculate_distance_from_five_cap(exons, strand, position):
    """Calculates distance from the 5' cap in a list of exons."""
    distance = 0
    sorted_exons = sorted(exons, key=lambda x: x[0])
    if strand == '-':
        sorted_exons = list(reversed(sorted_exons))
    for exon_start, exon_end in sorted_exons:
        if (strand == '+' and exon_end < position) or (strand == '-' and exon_start > position):
            distance += exon_end - exon_start + 1
        elif (strand == '+' and exon_start <= position <= exon_end) or (strand == '-' and exon_end >= position >= exon_start):
            distance += position - exon_start if strand == '+' else exon_end - position
            break
    return distance

def process_variant_spliceai_2(variant, UTRs_by_gene, Introns_by_transcript, cutoff):
    """Processes a single variant and returns the processed result."""
    CHR = variant[0] if 'chr' in variant[0] else 'chr' + variant[0]
    POS = int(variant[1])
    REF = variant[3]
    ALT = '' if variant[4] in ('<DEL>', '.') else variant[4]
    SPLICEAI = variant[-1]
    try:
        GENE, AG_score, AL_score, DG_score, DL_score, AG_POS, AL_POS, DG_POS, DL_POS = SPLICEAI.split('|')
        AG_POS, AL_POS, DG_POS, DL_POS = int(AG_POS) + POS, int(AL_POS) + POS, int(DG_POS) + POS, int(DL_POS) + POS
        AG_score, AL_score, DG_score, DL_score = map(float, (AG_score, AL_score, DG_score, DL_score))
    except ValueError:
        # If SPLICEAI field does not have expected format, skip processing
        return []
    gene_UTRs = UTRs_by_gene.get(GENE, [])
    result = []
    for UTR in gene_UTRs:
        if CHR != UTR[0] or int(UTR[1]) > POS or int(UTR[2]) < POS:
            continue
        exons = ast.literal_eval(UTR[13])
        strand = UTR[3]
        if strand == '+':
            # Process AG_score
            if AG_score >= cutoff:
                if all(AL_POS not in range(t[0], t[1]+1) for t in exons):
                    continue
                if AG_POS < AL_POS:
                    variant_type = 'AG_insertion_+'
                    transcript_INTRONs = Introns_by_transcript.get(UTR[6], [])
                    for intron in transcript_INTRONs:
                        if AL_POS == int(intron[2]):
                            newPOS = int(intron[1])
                            newREF = intron[11][0]
                            newALT = intron[11][0] + intron[11][AG_POS - AL_POS -1 : -1]
                            if AG_POS <= POS < AL_POS and AG_POS < POS + len(REF) -1 <= AL_POS:
                                newALT = newALT[:POS - AG_POS +1] + ALT + newALT[POS - AG_POS + len(REF) +1 :]
                            result.append([CHR, newPOS, variant[2], newREF, newALT] + variant[5:] +
                                          [UTR[6], f'{CHR}_{POS}_{variant[2]}_{REF}_{ALT}', variant_type])
                elif AG_POS > AL_POS:
                    variant_type = 'AG_deletion_+'
                    new = calculate_distance_from_five_cap(exons, strand, AG_POS)
                    old = calculate_distance_from_five_cap(exons, strand, AL_POS)
                    newPOS = next((exons[i-1][1] for i in range(1, len(exons)) if exons[i][0] == AL_POS), None)
                    if newPOS:
                        newREF = UTR[12][old -1 : new]
                        newALT = UTR[12][old -1]
                        result.append([CHR, newPOS, variant[2], newREF, newALT] + variant[5:] +
                                      [UTR[6], f'{CHR}_{POS}_{variant[2]}_{REF}_{ALT}', variant_type])
            # Process DG_score
            if DG_score >= cutoff:
                if all(DL_POS not in range(t[0], t[1]+1) for t in exons):
                    continue
                if DG_POS > DL_POS:
                    variant_type = 'DG_insertion_+'
                    transcript_INTRONs = Introns_by_transcript.get(UTR[6], [])
                    for intron in transcript_INTRONs:
                        if DL_POS == int(intron[1]):
                            newPOS = int(intron[1])
                            newREF = intron[11][0]
                            newALT = intron[11][: DG_POS - DL_POS +1]
                            if DL_POS <= POS <= DG_POS and DL_POS <= POS + len(REF) -1 <= DG_POS:
                                newALT = newALT[: POS - DL_POS ] + ALT + newALT[POS - DL_POS + len(REF) :]
                            result.append([CHR, newPOS, variant[2], newREF, newALT] + variant[5:] +
                                          [UTR[6], f'{CHR}_{POS}_{variant[2]}_{REF}_{ALT}', variant_type])
                elif DG_POS < DL_POS:
                    variant_type = 'DG_deletion_+'
                    new = calculate_distance_from_five_cap(exons, strand, DG_POS)
                    old = calculate_distance_from_five_cap(exons, strand, DL_POS)
                    newPOS = DG_POS
                    newREF = UTR[12][new : old +1]
                    newALT = UTR[12][new]
                    result.append([CHR, newPOS, variant[2], newREF, newALT] + variant[5:] +
                                  [UTR[6], f'{CHR}_{POS}_{variant[2]}_{REF}_{ALT}', variant_type])
        else:  # strand == '-'
            if AG_score >= cutoff:
                if all(AL_POS not in range(t[0], t[1]+1) for t in exons): 
                    continue
                if AG_POS > AL_POS: 
                    variant_type = 'AG_insertion_-'
                    transcript_INTRONs = Introns_by_transcript.get(UTR[6], [])
                    for intron in transcript_INTRONs:
                        if AL_POS == int(intron[1]): 
                            newPOS = int(intron[1]) 
                            newREF = rev_seq(intron[11][-1]) 
                            newALT = rev_seq(intron[11][AL_POS - AG_POS -1: ])
                            if AL_POS <= POS <= AG_POS and AL_POS < POS + len(REF) - 1 <= AG_POS: 
                                newALT = newALT[: POS - AL_POS] + ALT + newALT[ POS - AG_POS + len(REF) -1 :]
                            result.append([CHR, newPOS, variant[2], newREF, newALT] + variant[5:] +
                                          [UTR[6], f'{CHR}_{POS}_{variant[2]}_{REF}_{ALT}', variant_type])
                elif AG_POS < AL_POS: 
                    variant_type = 'AG_deletion_-'
                    new = calculate_distance_from_five_cap(exons, strand, AG_POS) 
                    old = calculate_distance_from_five_cap(exons, strand, AL_POS) 
                    newPOS = AG_POS
                    newREF = rev_seq(UTR[12][old : new +1]) 
                    newALT = rev_seq(UTR[12][new])
                    result.append([CHR, newPOS, variant[2], newREF, newALT] + variant[5:] +
                                      [UTR[6], f'{CHR}_{POS}_{variant[2]}_{REF}_{ALT}', variant_type])
            if DG_score >= cutoff:
                if all(DL_POS not in range(t[0], t[1]+1) for t in exons): 
                    continue
                if DG_POS < DL_POS: 
                    variant_type = 'DG_insertion_-'
                    transcript_INTRONs = Introns_by_transcript.get(UTR[6], [])
                    for intron in transcript_INTRONs:
                        if DL_POS == int(intron[2]): 
                            newPOS = int(intron[1]) 
                            newREF = rev_seq(intron[11][-1])
                            newALT = newREF + rev_seq(intron[11][1: DL_POS - DG_POS + 1]) 
                            if DG_POS <= POS < DL_POS and DG_POS <= POS + len(REF) - 1 < DL_POS: 
                                newALT = newALT[:POS - DG_POS +1] + ALT + newALT[POS - DG_POS + len(REF) +1:]
                            result.append([CHR, newPOS, variant[2], newREF, newALT] + variant[5:] +
                                          [UTR[6], f'{CHR}_{POS}_{variant[2]}_{REF}_{ALT}', variant_type])
                elif DG_POS > DL_POS:
                    variant_type = 'DG_deletion_-'
                    new = calculate_distance_from_five_cap(exons, strand, DG_POS)
                    old = calculate_distance_from_five_cap(exons, strand, DL_POS)
                    newPOS = next((exons[i-1][1] for i in range(1, len(exons)) if exons[i][0] == DL_POS), None)
                    if newPOS:
                        newALT = rev_seq(UTR[12][old+1])
                        newREF = newALT + rev_seq(UTR[12][new +1: old +1])
                        result.append([CHR, newPOS, variant[2], newREF, newALT] + variant[5:] +
                                      [UTR[6], f'{CHR}_{POS}_{variant[2]}_{REF}_{ALT}', variant_type])
    return result

def process_variants_spliceai_2(input_file, output_file_path, data_dir, cutoff):
    """Processes all variants and writes the results to the output file."""
    UTRS_FILE_PATH = os.path.join(os.path.expanduser(data_dir), '5UTRs.tsv')
    INTRONS_FILE_PATH = os.path.join(os.path.expanduser(data_dir), 'Introns.tsv')
    UTRs = load_tsv_data(UTRS_FILE_PATH)
    Introns = load_tsv_data(INTRONS_FILE_PATH)
    variants = load_tsv_data(input_file)
    UTRs_by_gene = defaultdict(list)
    Introns_by_transcript = defaultdict(list)

    # Group UTRs and Introns
    for UTR in UTRs[1:]:
        GENE = UTR[5]
        UTRs_by_gene[GENE].append(UTR)
    UTRs_header = UTRs[0]

    for Intron in Introns[1:]:
        TRANSCRIPT = Intron[7]
        Introns_by_transcript[TRANSCRIPT].append(Intron)

    results = []
    for variant in variants[1:]:
        if ',' in variant[4]:
            continue
        processed_variant = process_variant_spliceai_2(variant, UTRs_by_gene, Introns_by_transcript, cutoff=cutoff)
        if processed_variant:
            results.extend(processed_variant)
    fields = variants[0] + [UTRs_header[6], 'True_variant', 'type']

    with open(output_file_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(fields)
        writer.writerows(results)

# Optional main function
def main():
    import argparse
    parser = argparse.ArgumentParser(description='Process variants with SpliceAI annotations.')
    parser.add_argument('input_file', type=str, help='Path to the input variants file.')
    parser.add_argument('output_file', type=str, help='Path to the output file.')
    parser.add_argument('--data-dir', type=str, default=os.path.join(os.path.expanduser("~"), ".5ULTRA", "data"), help='Path to the data directory.')
    parser.add_argument('--cutoff', type=float, default=0.2, help='Cutoff value for scores.')
    args = parser.parse_args()

    # Load variants
    process_variants_spliceai_2(args.input_file, args.output_file, data_dir=args.data_dir, cutoff=args.cutoff)

if __name__ == "__main__":
    main()
