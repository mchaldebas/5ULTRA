#!/usr/bin/env python3

import csv
import ast
from math import nan
from collections import defaultdict
import sys
import subprocess

VARIANTS_FILE_PATH = sys.argv[1]
UTR_FILE_PATH = '../data/5UTRs.tsv'
UORF_FILE_PATH = '../data/uORFs.tsv'
OUTPUT_FILE_PATH = sys.argv[2]
BASE_PAIRING = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N', '*': '*'}
KOZAK_STRENGTH = {'Weak': 0, 'Adequate': 1, 'Strong': 2, '': nan}
STOP_CODONS = {'TAA', 'TAG', 'TGA'}

# function to reverse negative strand sequence
def rev_seq(fwd_seq):
    return ''.join(BASE_PAIRING[nuc] for nuc in reversed(fwd_seq))

def is_valid_number(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

# function to go from genomic positions to mRNA positions
def calculate_distance_from_five_cap(exons, strand, position):
    exons.sort(key=lambda x: x[0])
    distance = 0
    exons = exons if strand == '+' else reversed(exons)
    for exon_start, exon_end in exons:
        if (strand == '+' and exon_end < position) or (strand == '-' and exon_start > position):
            distance += exon_end - exon_start + 1
        elif (strand == '+' and exon_start < position) or (strand == '-' and exon_end > position):
            distance += position - exon_start if strand == '+' else exon_end - position
            break
    return distance

def calculate_genomic_position_from_five_cap(exons, strand, distance):
    exons.sort(key=lambda x: x[0])
    exons = exons if strand == '+' else list(reversed(exons))
    remaining_distance = distance
    for exon_start, exon_end in exons:
        exon_length = exon_end - exon_start + 1
        if remaining_distance <= exon_length:
            # Found the position within this exon
            return exon_start + remaining_distance if strand == '+' else exon_end - remaining_distance
        else:
            # Subtract the exon length from the remaining distance and move to the next exon
            remaining_distance -= exon_length
    return exon_start + remaining_distance if strand == '+' else exon_end - remaining_distance

# read tsv file into a list of list
def load_tsv_data(file_path):
    data = []
    try:
        with open(file_path, 'r') as file:
            csv_reader = csv.reader(file, delimiter='\t')
            for row in csv_reader:
                data.append(row)
    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
    return data

def load_vcf_data(vcf_path):
    data = []
    try :
        with open(vcf_path, 'r') as file:
            for line in file:
                if line.startswith('##'):
                    continue  # Skip meta-information lines
                elif line.startswith('#CHROM'):
                    parts = line.strip().split('\t')
                    data.append(parts)  # Add header to matrix
                else:
                    # These are the data lines; extract relevant information
                    parts = line.strip().split('\t')
                    row = parts[:9]  # Get the first 9 columns of standard VCF fields
                    genotypes = [part.split(':')[0] for part in parts[9:]]  # Extract genotypes for each sample
                    data.append(row + genotypes)  # Add this row to the matrix
    except FileNotFoundError:
        print(f"File not found: {vcf_path}")
    except Exception as e:
        print(f"Error reading file {vcf_path}: {e}")
    return data

# read kozak sequence to return its strength
def calculate_kozak_strength(kozak_sequence):
    if not kozak_sequence:
        return ''
    conditions = [kozak_sequence[1] in ['A', 'G'],
                  kozak_sequence[-2] == 'G']
    if all(conditions):
        return 'Strong'
    elif any(conditions):
        return 'Adequate'
    else:
        return 'Weak'

def get_score(chrom, pos, file_path):
    try:
        result = subprocess.check_output(
            f"tabix {file_path} {chrom}:{pos}-{pos}", shell=True).decode('utf-8').strip()
        if result:
            return result.split('\t')[-1]
    except subprocess.CalledProcessError:
        return None

# function for annotation of created uORF (uStart gain)
def uStart_gain(relativePosition, mutatedSequence, startPOS, STRAND, exons, CHR):
    uORF_START = relativePosition -2
    # search for START position
    while mutatedSequence[uORF_START: uORF_START +3] != 'ATG': uORF_START += 1
    # scan for inframe STOP
    uORF_END = uORF_START +3
    while mutatedSequence[uORF_END : uORF_END +3] not in STOP_CODONS and uORF_END < len(mutatedSequence):
        uORF_END += 3
    # annotation
    uORF_END += 2
    uSTART_mSTART_DIST = startPOS - uORF_START
    uSTOP_CODON = mutatedSequence[uORF_END -2 : uORF_END +1]
    uORF_TYPE = 'Non-Overlapping' if uORF_END < startPOS else 'N-terminal extension' if (startPOS - uORF_START) %3 == 0 else 'Overlapping'
    uKOZAK = mutatedSequence[uORF_START -4 : uORF_START +5]
    uKOZAK_STRENGTH = calculate_kozak_strength(uKOZAK)
    uORF_LENGTH = uORF_END - uORF_START +1 if uORF_TYPE != 'N-terminal extension' else startPOS - uORF_START
    if uORF_TYPE != 'N-terminal extension':
        uORF_SEQ = '' #mutatedSequence[uORF_START:uORF_END]
        uORF_rank = ''
    else:
        uORF_SEQ = '' #mutatedSequence[uORF_START:startPOS]
        uORF_rank = ''
    uORF_START = calculate_genomic_position_from_five_cap(exons, STRAND, uORF_START)
    if STRAND == '+':
        pos1, pos2, pos3 = uORF_START, uORF_START + 1, uORF_START + 2
    else:
        pos1, pos2, pos3 = uORF_START, uORF_START - 1, uORF_START - 2
    # Get Phylop scores for the 3 nucleotides
    phyloP_scores = [
        get_score(CHR, pos1, "./data/PhyloP/5UTR.hg38.phyloP100way/{}.bed.gz".format(CHR)),
        get_score(CHR, pos2, "./data/PhyloP/5UTR.hg38.phyloP100way/{}.bed.gz".format(CHR)),
        get_score(CHR, pos3, "./data/PhyloP/5UTR.hg38.phyloP100way/{}.bed.gz".format(CHR))
    ]
    # Calculate mean Phylop score
    phyloP_scores = [float(score) for score in phyloP_scores if score and (score.replace('.', '', 1).replace('-', '', 1).isdigit())]
    mean_phylop = sum(phyloP_scores) / len(phyloP_scores) if phyloP_scores else "NA"
    # Get PhastCons scores for the 3 nucleotides
    phastCons_scores = [
        get_score(CHR, pos1, "./data/Phastcons/5UTR.hg38.phastCons100way/{}.bed.gz".format(CHR)),
        get_score(CHR, pos2, "./data/Phastcons/5UTR.hg38.phastCons100way/{}.bed.gz".format(CHR)),
        get_score(CHR, pos3, "./data/Phastcons/5UTR.hg38.phastCons100way/{}.bed.gz".format(CHR))
    ]
    # Calculate mean PhastCons score
    phastCons_scores = [float(score) for score in phastCons_scores if score and score.replace('.', '', 1).isdigit()]
    mean_phastcons = sum(phastCons_scores) / len(phastCons_scores) if phastCons_scores else "NA"
    uORF_END = calculate_genomic_position_from_five_cap(exons, STRAND, uORF_END)
    return [uORF_START, uORF_END, '000', uSTART_mSTART_DIST, 'ATG', uSTOP_CODON, uORF_TYPE, uKOZAK, uKOZAK_STRENGTH, uORF_LENGTH, uORF_LENGTH/3, uORF_SEQ, uORF_rank, mean_phylop, mean_phastcons]

def process_variant(variant, utrs_by_transcript, uorfs_by_transcript):
    if not variant[1].isdigit():
        print(f"Warning: Skipping variant with invalid position value: {variant}")
        return None  # or handle it as needed
    CHR = variant[0] if 'chr' in variant[0] else 'chr' + variant[0]
    POS = int(variant[1])
    REF = variant[3]
    ALT = variant[4] if variant[4] != '<DEL>' else ''
    ALT = ALT if ALT != '.' else ''
    transcript_utrs = utrs_by_transcript.get(variant[-3], [])
    result = []
    for UTR in transcript_utrs:
        CSQ = [[],[]]
        uORFAnnotations = []
        # check if variant is on the same chromosome and in the 5UTR bondaries
        if not (int(UTR[1]) < POS < int(UTR[2])):
            continue
        exons = ast.literal_eval(UTR[13])
        # retreive relative position, wild type and mutated sequences
        relativePosition = calculate_distance_from_five_cap(exons, UTR[3], POS) if UTR[3] == '+' and len(REF) == 1 \
            else calculate_distance_from_five_cap(exons, UTR[3], POS +len(REF) -1)
        wtSEQ = UTR[12]
        mutatedSequence = wtSEQ[:relativePosition] + ALT + wtSEQ[relativePosition+len(REF):] if UTR[3] == '+' \
            else wtSEQ[:relativePosition] + rev_seq(ALT) + wtSEQ[relativePosition+len(REF):]
        startPOS = calculate_distance_from_five_cap(exons, UTR[3], int(UTR[2])) if UTR[3] == '+' \
            else calculate_distance_from_five_cap(exons, UTR[3], int(UTR[1]))
        startPOS += len(ALT) - len(REF)
        # exclude mStart loss variants
        if 'ATG' != mutatedSequence[startPOS :startPOS +3]:
            continue
        # mKozaks
        newKOZAK = mutatedSequence[startPOS -4 :startPOS +5]
        if len(newKOZAK) == 9 and len(UTR[10]) == 9: 
            if UTR[10][-2] != newKOZAK[-2] or UTR[10][1] != newKOZAK[1] or UTR[10][3] != newKOZAK[3]:
                # compare WT kozak strength with the Mutated kozak strength
                newKOZAK_STRENGTH = calculate_kozak_strength(newKOZAK)
                newKOZAK_STRENGTH = UTR[10] if newKOZAK_STRENGTH == '' else newKOZAK_STRENGTH
                if KOZAK_STRENGTH[newKOZAK_STRENGTH] < KOZAK_STRENGTH[UTR[11]]:
                    CSQ[0].extend(['mKozak'])
                    CSQ[1].extend(['decreased'])
                    uORFAnnotations += [['']*15]
                if KOZAK_STRENGTH[newKOZAK_STRENGTH] > KOZAK_STRENGTH[UTR[11]]:
                    CSQ[0].extend(['mKozak'])
                    CSQ[1].extend(['increased'])
                    uORFAnnotations += [['']*15]
        # uStart gain
        if 'ATG' in mutatedSequence[relativePosition-2: relativePosition+len(ALT)+2] and 'ATG' not in wtSEQ[ relativePosition-2: relativePosition+len(REF)+2]:
            CSQ[0].extend(['uStart_gain'])
            Anno = uStart_gain(relativePosition, mutatedSequence, startPOS, UTR[3], exons, CHR)
            uORFAnnotations += [Anno]
            if uORFAnnotations[-1][6] != 'N-terminal extension':
                CSQ[1].extend(['decreased'])
            else: CSQ[1].extend(['N-terminal extension'])
        # check if 5UTR has existing uORF(s)
        if float(UTR[14]) != 0:
            transcript_ids = UTR[6]
            transcript_uorfs = uorfs_by_transcript.get(transcript_ids,[])
            for uORF in transcript_uorfs:
                uSTART = int(uORF[8]) - int(uORF[17])
                uSTOP = uSTART + int(uORF[23]) -3
                if uSTOP >= relativePosition and len(REF) < len(ALT):
                    if uSTART >= relativePosition:
                        uSTART += len(ALT) -1
                    uSTOP += len(ALT) -1
                if uSTOP >= relativePosition +len(REF) and len(REF) > len(ALT):
                    if uSTART >= relativePosition +len(REF):
                        uSTART -= len(REF) -1
                    uSTOP -= len(REF) -1
                # check if variant is in uORF + uKozak
                if uSTART-6 <= relativePosition <= uSTOP+2:
                    # uStart loss & uKozak
                    if mutatedSequence[uSTART : uSTART +3] != uORF[18] and mutatedSequence[uSTART : uSTART +3] != 'ATG':
                        CSQ[0].extend(['uStart_loss'])
                        CSQ[1].extend(['increased'])
                        Anno = uORF[1:3] + [uORF[4]] + uORF[17:-4] + uORF[-3:]
                        uORFAnnotations += [Anno]
                        continue
                    # scan frame for STOP then uStop gain & uStop loss
                    codon = uSTART
                    while mutatedSequence[codon : codon +3] not in STOP_CODONS and codon < len(mutatedSequence):
                        codon += 3
                    NewUstopCodon = mutatedSequence[codon : codon+3]
                    if codon < uSTOP and codon +2 < startPOS:
                        if uORF[20] == 'Overlapping':
                            CSQ[0].extend(['uStop_gain to Non-Overlapping'])
                            CSQ[1].extend(['increased'])
                            Anno = uORF[1:3] + [uORF[4]] + uORF[17:19] + [uORF[19] + " > " + NewUstopCodon] + uORF[20:-4] + uORF[-3:]
                            uORFAnnotations += [Anno]
                            continue
                        elif uORF[20] == 'N-terminal extension':
                            CSQ[0].extend(['uStop_gain to Non-Overlapping'])
                            CSQ[1].extend(['decreased'])
                            Anno = uORF[1:3] + [uORF[4]] + uORF[17:19] + [uORF[19] + " > " + NewUstopCodon] + uORF[20:-4] + uORF[-3:]
                            uORFAnnotations += [Anno]
                        else:
                            CSQ[0].extend(['uStop_gain shorter Non-Overlapping'])
                            CSQ[1].extend(['increased'])
                            Anno = uORF[1:3] + [uORF[4]] + uORF[17:19] + [uORF[19] + " > " + NewUstopCodon] + uORF[20:-4] + uORF[-3:]
                            uORFAnnotations += [Anno]
                            continue
                    if codon > uSTOP and uORF[20] == 'Non-Overlapping':
                        if (codon - startPOS)%3 == 0:
                            CSQ[0].extend(['uStop_loss to N-terminal extension'])
                            CSQ[1].extend(['N-terminal extension'])
                            Anno = uORF[1:3] + [uORF[4]] + uORF[17:19] + [uORF[19] + " > " + NewUstopCodon] + uORF[20:-4] + uORF[-3:]
                            uORFAnnotations += [Anno]
                        elif codon +2 > startPOS:
                            CSQ[0].extend(['uStop_loss to Overlapping'])
                            CSQ[1].extend(['decreased'])
                            Anno = uORF[1:3] + [uORF[4]] + uORF[17:19] + [uORF[19] + " > " + NewUstopCodon] + uORF[20:-4] + uORF[-3:]
                            uORFAnnotations += [Anno]
                        else:
                            CSQ[0].extend(['uStop_loss longer Non-Overlapping'])
                            CSQ[1].extend(['decreased'])
                            Anno = uORF[1:3] + [uORF[4]] + uORF[17:19] + [uORF[19] + " > " + NewUstopCodon] + uORF[20:-4] + uORF[-3:]
                            uORFAnnotations += [Anno]
                    if uSTART -1 == relativePosition or relativePosition == uSTART -3 or relativePosition == uSTART +3:
                        newKOZAK = mutatedSequence[uSTART -4 :uSTART +5]
                        # compare WT kozak strength with the Mutated kozak strength
                        newKOZAK_STRENGTH = calculate_kozak_strength(newKOZAK)
                        newKOZAK_STRENGTH = uORF[22] if newKOZAK_STRENGTH == '' else newKOZAK_STRENGTH
                        if KOZAK_STRENGTH[newKOZAK_STRENGTH] < KOZAK_STRENGTH[uORF[22]]:
                            CSQ[0].extend(['uKozak'])
                            CSQ[1].extend(['increased'])
                            Anno = uORF[1:3] + [uORF[4]] + uORF[17:-4] + uORF[-3:]
                            uORFAnnotations += [Anno]
                        elif KOZAK_STRENGTH[newKOZAK_STRENGTH] > KOZAK_STRENGTH[uORF[22]]:
                            CSQ[0].extend(['uKozak'])
                            CSQ[1].extend(['decreased'])
                            Anno= uORF[1:3] + [uORF[4]] + uORF[17:-4] + uORF[-3:]
                            uORFAnnotations += [Anno]
        count = 0
        for hit in CSQ[0]:
            result.append( variant[-2].split('_') + variant[5:-4] + [variant[-4], variant[-1]] + [hit, CSQ[1][count]] + UTR[1:12] + UTR[14:] + uORFAnnotations[count])
            count += 1
    return result

def main():
    if VARIANTS_FILE_PATH[-4:] == '.vcf':
        variants = load_vcf_data(VARIANTS_FILE_PATH)
    else:
        variants = load_tsv_data(VARIANTS_FILE_PATH)
    UTRs = load_tsv_data(UTR_FILE_PATH)
    uORFs = load_tsv_data(UORF_FILE_PATH)
    utrs_by_transcript = defaultdict(list)
    uorfs_by_transcript = defaultdict(list)
    for UTR in UTRs[1:]:
        TRANSCRIPTS = UTR[6]
        utrs_by_transcript[TRANSCRIPTS].append(UTR)
    UTRs = UTRs[0]
    for uORF in uORFs[1:]:
        TRANSCRIPTS = uORF[5]
        uorfs_by_transcript[TRANSCRIPTS].append(uORF)
    uORFs = uORFs[0]
    with open(OUTPUT_FILE_PATH, 'w') as f:
        write = csv.writer(f, delimiter='\t')
        fields = variants[0][:-4] + ['SpliceAI', 'variant_type', 'CSQ', 'translation'] + UTRs[1:12] + UTRs[14:] + uORFs[1:3] + [uORFs[4]] + uORFs[17:-4] + uORFs[-3:]
        write.writerow(fields)
        for variant in variants[1:]:
            if ',' in variant[4]:
                continue
            processed_variant = process_variant(variant, utrs_by_transcript, uorfs_by_transcript)
            if processed_variant:
                write.writerows(processed_variant)

if __name__ == "__main__":
    main()
