import csv
import ast
from math import nan
from collections import defaultdict
import pysam
import os

BASE_PAIRING = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N', '*': '*'}
KOZAK_STRENGTH = {'Weak': 0, 'Adequate': 1, 'Strong': 2, '': nan}
STOP_CODONS = {'TAA', 'TAG', 'TGA'}

def rev_seq(fwd_seq):
    """Reverse complements a DNA sequence."""
    return ''.join(BASE_PAIRING.get(nuc, 'N') for nuc in reversed(fwd_seq))

def is_valid_number(value):
    """Checks if a value can be converted to a float."""
    try:
        float(value)
        return True
    except ValueError:
        return False

def calculate_distance_from_five_cap(exons, strand, position):
    """Calculates the distance from the 5' cap to a given genomic position."""
    exons.sort(key=lambda x: x[0])
    distance = 0
    exons = exons if strand == '+' else list(reversed(exons))
    for exon_start, exon_end in exons:
        if (strand == '+' and exon_end < position) or (strand == '-' and exon_start > position):
            distance += exon_end - exon_start + 1
        elif (strand == '+' and exon_start < position) or (strand == '-' and exon_end > position):
            distance += position - exon_start if strand == '+' else exon_end - position
            break
    return distance

def calculate_genomic_position_from_five_cap(exons, strand, distance):
    """Calculates the genomic position from a given distance from the 5' cap."""
    exons.sort(key=lambda x: x[0])
    exons = exons if strand == '+' else list(reversed(exons))
    remaining_distance = distance
    for exon_start, exon_end in exons:
        exon_length = exon_end - exon_start + 1
        if remaining_distance <= exon_length:
            return exon_start + remaining_distance if strand == '+' else exon_end - remaining_distance
        else:
            remaining_distance -= exon_length
    return exon_start + remaining_distance if strand == '+' else exon_end - remaining_distance

def load_tsv_data(file_path):
    """Loads TSV data from a file."""
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
    """Loads VCF data from a file."""
    data = []
    try:
        with open(vcf_path, 'r') as file:
            for line in file:
                if line.startswith('##'):
                    continue
                elif line.startswith('#CHROM'):
                    parts = line.strip().split('\t')
                    data.append(parts)
                else:
                    parts = line.strip().split('\t')
                    row = parts[:9]
                    genotypes = [part.split(':')[0] for part in parts[9:]]
                    data.append(row + genotypes)
    except FileNotFoundError:
        print(f"File not found: {vcf_path}")
    except Exception as e:
        print(f"Error reading file {vcf_path}: {e}")
    return data

def calculate_kozak_strength(kozak_sequence):
    """Determines the Kozak sequence strength."""
    if not kozak_sequence or len(kozak_sequence) < 9:
        return ''
    conditions = [kozak_sequence[1] in ['A', 'G'], kozak_sequence[-2] == 'G']
    if all(conditions):
        return 'Strong'
    elif any(conditions):
        return 'Adequate'
    else:
        return 'Weak'

def get_score(chrom, pos, file_path):
    """Fetches conservation scores from tabix-indexed files."""
    try:
        tabix_file = pysam.TabixFile(file_path)
        records = tabix_file.fetch(chrom, pos - 1, pos)
        for record in records:
            return record.split('\t')[-1]
        return None
    except (OSError, ValueError, KeyError, pysam.TabixError):
        return None

def uStart_gain(relativePosition, mutatedSequence, startPOS, STRAND, exons, CHR):
    """Annotates created uORFs (uStart gain)."""
    uORF_START = relativePosition - 2
    while mutatedSequence[uORF_START: uORF_START + 3] != 'ATG':
        uORF_START += 1
    uORF_END = uORF_START + 3
    while mutatedSequence[uORF_END: uORF_END + 3] not in STOP_CODONS and uORF_END < len(mutatedSequence):
        uORF_END += 3
    uORF_END += 2
    uSTART_mSTART_DIST = startPOS - uORF_START
    uSTOP_CODON = mutatedSequence[uORF_END - 2: uORF_END + 1]
    uORF_TYPE = 'Non-Overlapping' if uORF_END < startPOS else 'N-terminal extension' if (startPOS - uORF_START) % 3 == 0 else 'Overlapping'
    uKOZAK = mutatedSequence[uORF_START - 4: uORF_START + 5]
    uKOZAK_STRENGTH = calculate_kozak_strength(uKOZAK)
    uORF_LENGTH = uORF_END - uORF_START + 1 if uORF_TYPE != 'N-terminal extension' else startPOS - uORF_START
    uORF_START_GENOMIC = calculate_genomic_position_from_five_cap(exons, STRAND, uORF_START)
    if STRAND == '+':
        pos1, pos2, pos3 = uORF_START_GENOMIC, uORF_START_GENOMIC + 1, uORF_START_GENOMIC + 2
    else:
        pos1, pos2, pos3 = uORF_START_GENOMIC, uORF_START_GENOMIC - 1, uORF_START_GENOMIC - 2
    phyloP_scores = [
        get_score(CHR, pos1, f"./data/5UTR.hg38.phyloP100way/{CHR}.bed.gz"),
        get_score(CHR, pos2, f"./data/5UTR.hg38.phyloP100way/{CHR}.bed.gz"),
        get_score(CHR, pos3, f"./data/5UTR.hg38.phyloP100way/{CHR}.bed.gz")
    ]
    phyloP_scores = [float(score) for score in phyloP_scores if score and is_valid_number(score)]
    mean_phylop = sum(phyloP_scores) / len(phyloP_scores) if phyloP_scores else "NA"
    phastCons_scores = [
        get_score(CHR, pos1, f"./data/5UTR.hg38.phastCons100way/{CHR}.bed.gz"),
        get_score(CHR, pos2, f"./data/5UTR.hg38.phastCons100way/{CHR}.bed.gz"),
        get_score(CHR, pos3, f"./data/5UTR.hg38.phastCons100way/{CHR}.bed.gz")
    ]
    phastCons_scores = [float(score) for score in phastCons_scores if score and is_valid_number(score)]
    mean_phastcons = sum(phastCons_scores) / len(phastCons_scores) if phastCons_scores else "NA"
    uORF_END_GENOMIC = calculate_genomic_position_from_five_cap(exons, STRAND, uORF_END)
    return [
        uORF_START_GENOMIC, uORF_END_GENOMIC, '000', uSTART_mSTART_DIST, 'ATG',
        uSTOP_CODON, uORF_TYPE, uKOZAK, uKOZAK_STRENGTH,
        uORF_LENGTH, uORF_LENGTH / 3, '', '', mean_phylop, mean_phastcons
    ]

def process_variant(variant, utrs_by_chromosome, uorfs_by_transcript):
    """Processes a single variant and returns the annotated result."""
    CHR = variant[0] if 'chr' in variant[0] else 'chr' + variant[0]
    POS = int(variant[1])
    REF = variant[3]
    ALT = variant[4] if variant[4] != '<DEL>' else ''
    ALT = ALT if ALT != '.' else ''
    chromosome_utrs = utrs_by_chromosome.get(CHR, [])
    result = []
    for UTR in chromosome_utrs:
        CSQ = [[], []]
        uORFAnnotations = []
        if UTR[0] != CHR or not (int(UTR[1]) < POS < int(UTR[2])):
            continue
        status = 'out'
        exons = ast.literal_eval(UTR[13])
        for exon in exons:
            if exon[0] <= POS <= exon[1] and exon[0] <= POS + len(REF) - 1 <= exon[1]:
                status = 'in'
                break
        if status == 'out':
            continue
        relativePosition = calculate_distance_from_five_cap(exons, UTR[3], POS) if UTR[3] == '+' and len(REF) == 1 \
            else calculate_distance_from_five_cap(exons, UTR[3], POS + len(REF) - 1)
        wtSEQ = UTR[12]
        mutatedSequence = wtSEQ[:relativePosition] + (ALT if UTR[3] == '+' else rev_seq(ALT)) + wtSEQ[relativePosition+len(REF):]
        startPOS = calculate_distance_from_five_cap(exons, UTR[3], int(UTR[2])) if UTR[3] == '+' \
            else calculate_distance_from_five_cap(exons, UTR[3], int(UTR[1]))
        startPOS += len(ALT) - len(REF)
        if 'ATG' != mutatedSequence[startPOS:startPOS+3]:
            continue
        newKOZAK = mutatedSequence[startPOS - 4:startPOS + 5]
        if len(newKOZAK) == 9 and len(UTR[10]) == 9:
            if UTR[10][-2] != newKOZAK[-2] or UTR[10][1] != newKOZAK[1] or UTR[10][3] != newKOZAK[3]:
                newKOZAK_STRENGTH = calculate_kozak_strength(newKOZAK)
                newKOZAK_STRENGTH = UTR[10] if newKOZAK_STRENGTH == '' else newKOZAK_STRENGTH
                if KOZAK_STRENGTH[newKOZAK_STRENGTH] < KOZAK_STRENGTH[UTR[11]]:
                    CSQ[0].append('mKozak')
                    CSQ[1].append('decreased')
                    uORFAnnotations.append([''] * 15)
                if KOZAK_STRENGTH[newKOZAK_STRENGTH] > KOZAK_STRENGTH[UTR[11]]:
                    CSQ[0].append('mKozak')
                    CSQ[1].append('increased')
                    uORFAnnotations.append([''] * 15)
        if 'ATG' in mutatedSequence[relativePosition-2: relativePosition+len(ALT)+2] and 'ATG' not in wtSEQ[relativePosition-2: relativePosition+len(REF)+2]:
            CSQ[0].append('uStart_gain')
            Anno = uStart_gain(relativePosition, mutatedSequence, startPOS, UTR[3], exons, CHR)
            uORFAnnotations.append(Anno)
            CSQ[1].append('N-terminal extension' if uORFAnnotations[-1][6] == 'N-terminal extension' else 'decreased')
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
                        if codon >= startPOS and (codon - startPOS) %3 == 0:
                            CSQ[0].extend(['uStop_loss to N-terminal extension'])
                            CSQ[1].extend(['N-terminal extension'])
                            Anno = uORF[1:3] + [uORF[4]] + uORF[17:19] + [uORF[19] + " > " + NewUstopCodon] + uORF[20:-4] + uORF[-3:]
                            uORFAnnotations += [Anno]
                        elif codon >= startPOS:
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
            result.append(variant + [hit, CSQ[1][count]] + UTR[1:12] + UTR[14:] + uORFAnnotations[count])
            count += 1
    return result

def process_variants(input_variants, output_file_path, data_dir='./data'):
    """Processes all variants and writes the results to the output file."""
    UTR_FILE_PATH = os.path.join(data_dir, '5UTRs.tsv')
    UORF_FILE_PATH = os.path.join(data_dir, 'uORFs.tsv')
    UTRs = load_tsv_data(UTR_FILE_PATH)
    uORFs = load_tsv_data(UORF_FILE_PATH)
    utrs_by_chromosome = defaultdict(list)
    uorfs_by_transcript = defaultdict(list)
    for UTR in UTRs[1:]:
        CHR = UTR[0]
        utrs_by_chromosome[CHR].append(UTR)
    UTR_headers = UTRs[0]
    for uORF in uORFs[1:]:
        TRANSCRIPTS = uORF[5]
        uorfs_by_transcript[TRANSCRIPTS].append(uORF)
    uORF_headers = uORFs[0]
    with open(output_file_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        fields = input_variants[0] + ['CSQ', 'translation'] + UTR_headers[1:12] + UTR_headers[14:] + uORF_headers[1:3] + [uORF_headers[4]] + uORF_headers[17:-4] + uORF_headers[-3:]
        writer.writerow(fields)
        for variant in input_variants[1:]:
            if ',' in variant[4]:
                continue
            processed_variant = process_variant(variant, utrs_by_chromosome, uorfs_by_transcript)
            if processed_variant:
                writer.writerows(processed_variant)

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Detect variants.')
    parser.add_argument('input_file_path', type=str, help='Path to the filtered input file.')
    parser.add_argument('output_file_path', type=str, help='Path to the detection output file.')
    parser.add_argument('--data-dir', type=str, default='./data', help='Path to the data directory.')
    args = parser.parse_args()
    if args.input_file_path.endswith('.vcf'):
        variants = load_vcf_data(args.input_file_path)
    else:
        variants = load_tsv_data(args.input_file_path)
    process_variants(variants, args.output_file_path, data_dir=args.data_dir)

if __name__ == "__main__":
    main()
