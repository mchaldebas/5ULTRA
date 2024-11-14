#!/usr/bin/env python3

import sys
import csv
from collections import defaultdict
import ast

# File paths
VARIANTS_FILE_PATH = sys.argv[1]
UTRS_FILE_PATH = './data/5UTRs.tsv'
INTRONS_FILE_PATH = './data/52.filtered-Introns.tsv'
OUTPUT_FILE_PATH = sys.argv[2]

# Constants
BASE_PAIRING = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N', '*': '*'}
CUTOFF = 0.2

# Reverse sequence based on BASE_PAIRING dictionary
def rev_seq(fwd_seq):
    return ''.join(BASE_PAIRING[nuc] for nuc in reversed(fwd_seq))

# Load TSV data from a file
def load_tsv_data(file_path):
    data = []
    try:
        max_int = sys.maxsize
        while True:
            try:
                csv.field_size_limit(max_int)
                break
            except OverflowError:
                max_int //= 10
        with open(file_path, 'r') as file:
            csv_reader = csv.reader(file, delimiter='\t')
            data = list(csv_reader)
    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
    return data

# Calculate distance from the five cap in a list of exons
def calculate_distance_from_five_cap(exons, strand, position):
    distance = 0
    sorted_exons = sorted(exons, key=lambda x: x[0])
    if strand == '-':
        sorted_exons = list(reversed(sorted_exons))
    for exon_start, exon_end in sorted_exons:
        if (strand == '+' and exon_end < position) or (strand == '-' and exon_start > position):
            distance += exon_end - exon_start + 1
        elif (strand == '+' and exon_start < position) or (strand == '-' and exon_end > position):
            distance += position - exon_start if strand == '+' else exon_end - position
            break
    return distance

def process_variant(variant, UTRs_by_gene, Introns_by_transcript):
    CHR = variant[0] if 'chr' in variant[0] else 'chr' + variant[0]
    POS = int(variant[1])
    REF = variant[3]
    ALT = '' if variant[4] in ('<DEL>', '.') else variant[4]
    SPLICEAI = variant[-1]
    GENE, AG_score, AL_score, DG_score, DL_score, AG_POS, AL_POS, DG_POS, DL_POS = SPLICEAI.split('|')
    AG_POS, AL_POS, DG_POS, DL_POS = int(AG_POS) +POS, int(AL_POS) +POS, int(DG_POS) +POS, int(DL_POS) +POS
    AG_score, AL_score, DG_score, DL_score = map(float, (AG_score, AL_score, DG_score, DL_score))
    gene_UTRs = UTRs_by_gene.get(GENE, [])
    result = []
    for UTR in gene_UTRs:
        if CHR != UTR[0] or int(UTR[1]) > POS or int(UTR[2]) < POS:
            continue
        exons = ast.literal_eval(UTR[13])
        strand = UTR[3]
        if strand == '+':
            if AG_score >= CUTOFF:
                if all(AL_POS not in t for t in exons):
                    continue
                if AG_POS < AL_POS:
                    variant_type = 'AG_insertion_+'
                    transcript_INTRONs = Introns_by_transcript.get(UTR[6], [])
                    for intron in transcript_INTRONs:
                        if AL_POS == int(intron[2]):
                            newPOS = int(intron[1])
                            newREF = intron[11][0]
                            newALT = intron[11][0] + intron[11][AG_POS -AL_POS -1: -1]
                            if AG_POS <= POS <= AL_POS and AG_POS <= POS +len(REF)-1 <= AL_POS:
                                newALT = newALT[:POS -AG_POS] + ALT + newALT[POS -AG_POS + len(REF)-1:]
                            result.append([CHR, newPOS, variant[2], newREF, newALT] + variant[5:] + [UTR[6], f'{CHR}_{POS}_{variant[2]}_{REF}_{ALT}', variant_type])
                elif AG_POS > AL_POS:
                    variant_type = 'AG_deletion_+'
                    new = calculate_distance_from_five_cap(exons, strand, AG_POS)
                    old = calculate_distance_from_five_cap(exons, strand, AL_POS)
                    newPOS = next((exons[i-1][1] for i in range(1, len(exons)) if exons[i][0] == AL_POS), None)
                    if newPOS:
                        newREF = UTR[12][old -1: new]
                        newALT = UTR[12][old -1]
                        result.append([CHR, newPOS, variant[2], newREF, newALT] + variant[5:] + [UTR[6], f'{CHR}_{POS}_{variant[2]}_{REF}_{ALT}', variant_type])
            if DG_score >= CUTOFF:
                if all(DL_POS not in t for t in exons):
                    continue
                if DG_POS > DL_POS:
                    variant_type = 'DG_insertion_+'
                    transcript_INTRONs = Introns_by_transcript.get(UTR[6], [])
                    for intron in transcript_INTRONs:
                        if DL_POS == int(intron[1]):
                            newPOS = int(intron[1])
                            newREF = intron[11][0]
                            newALT = intron[11][:DG_POS -DL_POS +1]
                            if DL_POS <= POS <= DG_POS and DL_POS <= POS +len(REF)-1 <= DG_POS:
                                newALT = newALT[:POS -DL_POS] + ALT + newALT[POS -DL_POS + len(REF)-1:]
                            result.append([CHR, newPOS, variant[2], newREF, newALT] + variant[5:] + [UTR[6], f'{CHR}_{POS}_{variant[2]}_{REF}_{ALT}', variant_type])
                elif DG_POS < DL_POS:
                    variant_type = 'DG_deletion_+'
                    new = calculate_distance_from_five_cap(exons, strand, DG_POS)
                    old = calculate_distance_from_five_cap(exons, strand, DL_POS)
                    newPOS = DG_POS
                    newREF = UTR[12][new: old +1]
                    newALT = UTR[12][new]
                    result.append([CHR, newPOS, variant[2], newREF, newALT] + variant[5:] + [UTR[6], f'{CHR}_{POS}_{variant[2]}_{REF}_{ALT}', variant_type])
        else: # strand -
            if AG_score >= CUTOFF:
                if all(AL_POS not in t for t in exons):
                    continue
                if AG_POS > AL_POS:
                    variant_type = 'AG_insertion_-'
                    transcript_INTRONs = Introns_by_transcript.get(UTR[6], [])
                    for intron in transcript_INTRONs:
                        if AL_POS == int(intron[1]):
                            newPOS = int(intron[1])
                            newREF = rev_seq(intron[11][-1])
                            newALT = rev_seq(intron[11][AL_POS -AG_POS -1:])
                            if AG_POS >= POS >= AL_POS and AG_POS >= POS +len(REF)-1 >= AL_POS:
                                newALT = newALT[:POS -AL_POS] + ALT + newALT[POS -AL_POS + len(REF)-1:]
                            result.append([CHR, newPOS, variant[2], newREF, newALT] + variant[5:] + [UTR[6], f'{CHR}_{POS}_{variant[2]}_{REF}_{ALT}', variant_type])
                elif AG_POS < AL_POS:
                    variant_type = 'AG_deletion_-'
                    new = calculate_distance_from_five_cap(exons, strand, AG_POS)
                    old = calculate_distance_from_five_cap(exons, strand, AL_POS)
                    newPOS = AG_POS
                    newREF = rev_seq(UTR[12][old: new +1])
                    newALT = rev_seq(UTR[12][new])
                    result.append([CHR, newPOS, variant[2], newREF, newALT] + variant[5:] + [UTR[6], f'{CHR}_{POS}_{variant[2]}_{REF}_{ALT}', variant_type])
            if DG_score >= CUTOFF:
                if all(DL_POS not in t for t in exons):
                    continue
                if DG_POS < DL_POS:
                    variant_type = 'DG_insertion_-'
                    transcript_INTRONs = Introns_by_transcript.get(UTR[6], [])
                    for intron in transcript_INTRONs:
                        if DL_POS == int(intron[2]):
                            newPOS = int(intron[1])
                            newREF = rev_seq(intron[11][-1])
                            newALT = newREF + rev_seq(intron[11][1:DL_POS -DG_POS +1])
                            if DL_POS >= POS >= DG_POS and DL_POS >= POS +len(REF)-1 >= DG_POS:
                                newALT = newALT[:POS -DG_POS +1] + ALT + newALT[POS -DG_POS + len(REF) +1:]
                            result.append([CHR, newPOS, variant[2], newREF, newALT] + variant[5:] + [UTR[6], f'{CHR}_{POS}_{variant[2]}_{REF}_{ALT}', variant_type])
                elif DG_POS > DL_POS:
                    variant_type = 'DG_deletion_-'
                    new = calculate_distance_from_five_cap(exons, strand, DG_POS)
                    old = calculate_distance_from_five_cap(exons, strand, DL_POS)
                    newPOS = next((exons[i-1][1] for i in range(1, len(exons)) if exons[i][0] == DL_POS), None)
                    newREF = rev_seq(UTR[12][new +1: old +2])
                    newALT = rev_seq(UTR[12][old +1])
                    result.append([CHR, newPOS, variant[2], newREF, newALT] + variant[5:] + [UTR[6], f'{CHR}_{POS}_{variant[2]}_{REF}_{ALT}', variant_type])
    return result

def main():
    variants = load_tsv_data(VARIANTS_FILE_PATH)
    UTRs = load_tsv_data(UTRS_FILE_PATH)
    Introns = load_tsv_data(INTRONS_FILE_PATH)
    UTRs_by_gene = defaultdict(list)
    Introns_by_transcript = defaultdict(list)
    # Group UTRs and Introns by chromosome
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
        processed_variant = process_variant(variant, UTRs_by_gene, Introns_by_transcript)
        if processed_variant:
            results.extend(processed_variant)
    fields = variants[0] + [UTRs_header[6], 'True_variant', 'type']
    with open(OUTPUT_FILE_PATH, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(fields)
        writer.writerows(results)

if __name__ == "__main__":
    main()
