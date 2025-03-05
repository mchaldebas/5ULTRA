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
    try:
        with pysam.TabixFile(file_path) as tabix_file:  
            records = tabix_file.fetch(chrom, pos - 1, pos)
            for record in records:
                return record.split('\t')[-1]
            return None
    except OSError as e:
        print(f"OSError accessing file: {e}")
        return None
    except ValueError as e:
        print(f"ValueError (likely malformed input): {e}")
        return None
    except KeyError as e:
        print(f"KeyError (likely missing chromosome): {e}")
        return None
    except pysam.utils.SamtoolsError as e: 
        print(f"Pysam/Samtools error: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error: {e}")
        return None

def uStart_gain(relativePosition, mutatedSequence, startPOS, STRAND, exons, CHR, data_dir, POS, type, wtSEQ):
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
    if type == 'DG_insertion_+':
        uORF_START_GENOMIC = POS + (uORF_START - calculate_distance_from_five_cap(exons, STRAND, POS))
        uORF_END_GENOMIC = calculate_genomic_position_from_five_cap(exons, STRAND, uORF_END)
    elif type == 'DG_insertion_-':
        POS_exon = 0
        while exons[POS_exon][1] < POS:
            POS_exon += 1
        POS_exon += 1
        intron_end = exons[POS_exon][0]
        uORF_START_GENOMIC = POS + (intron_end - POS - uORF_START+ calculate_distance_from_five_cap(exons, STRAND, POS)) -1
        uORF_END_GENOMIC = calculate_genomic_position_from_five_cap(exons, STRAND, uORF_END)
    elif type == 'AG_insertion_+':
        POS_exon = 0
        while exons[POS_exon][1] < POS:
            POS_exon += 1
        POS_exon += 1
        intron_end = exons[POS_exon][0]
        uORF_START_GENOMIC = POS + (intron_end -POS +uORF_START -calculate_distance_from_five_cap(exons, STRAND, POS) - (len(mutatedSequence) - len(wtSEQ))) -1
        uORF_END_GENOMIC = calculate_genomic_position_from_five_cap(exons, STRAND, uORF_END)
    elif type == 'AG_insertion_-':
        uORF_START_GENOMIC = POS + (len(mutatedSequence) - len(wtSEQ) -uORF_START +calculate_distance_from_five_cap(exons, STRAND, POS))
        uORF_END_GENOMIC = calculate_genomic_position_from_five_cap(exons, STRAND, uORF_END)
    else: 
        uORF_START_GENOMIC = calculate_genomic_position_from_five_cap(exons, STRAND, uORF_START)
        uORF_END_GENOMIC = calculate_genomic_position_from_five_cap(exons, STRAND, uORF_END)
    if STRAND == '+':
        pos1, pos2, pos3 = uORF_START_GENOMIC, uORF_START_GENOMIC + 1, uORF_START_GENOMIC + 2
    else:
        pos1, pos2, pos3 = uORF_START_GENOMIC, uORF_START_GENOMIC - 1, uORF_START_GENOMIC - 2
    phyloP_scores = [
        get_score(CHR, pos1, os.path.join(os.path.expanduser(data_dir), '5UTR.hg38.phyloP100way', f"{CHR}.bed.gz")),
        get_score(CHR, pos2, os.path.join(os.path.expanduser(data_dir), '5UTR.hg38.phyloP100way', f"{CHR}.bed.gz")),
        get_score(CHR, pos3, os.path.join(os.path.expanduser(data_dir), '5UTR.hg38.phyloP100way', f"{CHR}.bed.gz"))
    ]
    phyloP_scores = [float(score) for score in phyloP_scores if score and is_valid_number(score)]
    mean_phylop = sum(phyloP_scores) / len(phyloP_scores) if phyloP_scores else "NA"
    phastCons_scores = [
        get_score(CHR, pos1, os.path.join(os.path.expanduser(data_dir), '5UTR.hg38.phastCons100way', f"{CHR}.bed.gz")),
        get_score(CHR, pos2, os.path.join(os.path.expanduser(data_dir), '5UTR.hg38.phastCons100way', f"{CHR}.bed.gz")),
        get_score(CHR, pos3, os.path.join(os.path.expanduser(data_dir), '5UTR.hg38.phastCons100way', f"{CHR}.bed.gz"))
    ]
    phastCons_scores = [float(score) for score in phastCons_scores if score and is_valid_number(score)]
    mean_phastcons = sum(phastCons_scores) / len(phastCons_scores) if phastCons_scores else "NA"
    return [
        uORF_START_GENOMIC, uORF_END_GENOMIC, '000', uSTART_mSTART_DIST, 'ATG',
        uSTOP_CODON, uORF_TYPE, uKOZAK, uKOZAK_STRENGTH,
        uORF_LENGTH, uORF_LENGTH / 3, '', '', mean_phylop, mean_phastcons
    ]

def process_variant_spliceai_3(variant, utrs_by_transcript, uorfs_by_transcript, data_dir):
    if not variant[1].isdigit():
        print(f"Warning: Skipping variant with invalid position value: {variant}")
        return None
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
        # check if variant is in the 5UTR bondaries
        if not (int(UTR[1]) <= POS <= int(UTR[2])):
            continue
        exons = ast.literal_eval(UTR[13])
        # retreive relative position, wild type and mutated sequences
        if UTR[3] == '+':
            relativePosition = calculate_distance_from_five_cap(exons, UTR[3], POS)
        else:
            relativePosition = calculate_distance_from_five_cap(exons, UTR[3], POS + len(REF) - 1)
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
            Anno = uStart_gain(relativePosition, mutatedSequence, startPOS, UTR[3], exons, CHR, data_dir, POS, variant[-1], wtSEQ)
            uORFAnnotations += [Anno]
            if uORFAnnotations[-1][6] != 'N-terminal extension':
                CSQ[1].extend(['decreased'])
            else: CSQ[1].extend(['N-terminal extension'])
        elif relativePosition < 2 and 'ATG' in mutatedSequence[: relativePosition+len(ALT)+2] and 'ATG' not in wtSEQ[: relativePosition+len(REF)+2]:
            CSQ[0].append('uStart_gain')
            Anno = uStart_gain(relativePosition, mutatedSequence, startPOS, UTR[3], exons, CHR, data_dir, POS, variant[-1], wtSEQ)
            uORFAnnotations.append(Anno)
            CSQ[1].append('N-terminal extension' if uORFAnnotations[-1][6] == 'N-terminal extension' else 'decreased')
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
                if (uSTART-6 <= relativePosition <= uSTOP+2) or (relativePosition <= uSTART <= relativePosition + len(REF)):
                    # uStart loss & uKozak
                    if mutatedSequence[uSTART : uSTART +3] != uORF[18] and mutatedSequence[uSTART : uSTART +3] != 'ATG':
                        CSQ[0].extend(['uStart_loss'])
                        CSQ[1].extend(['increased'])
                        Anno = uORF[1:3] + [uORF[4]] + uORF[17:-4] + uORF[-3:]
                        uORFAnnotations += [Anno]
                        continue
                    # scan frame for STOP then uStop gain & uStop loss
                    codon = uSTART
                    while mutatedSequence[codon : codon +3] not in STOP_CODONS and codon < len(mutatedSequence) and codon != startPOS:
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
                        elif uORF[20] == 'Non-Overlapping':
                            CSQ[0].extend(['uStop_gain shorter Non-Overlapping'])
                            CSQ[1].extend(['increased'])
                            Anno = uORF[1:3] + [uORF[4]] + uORF[17:19] + [uORF[19] + " > " + NewUstopCodon] + uORF[20:-4] + uORF[-3:]
                            uORFAnnotations += [Anno]
                            continue
                    elif codon < uSTOP and codon == startPOS and uORF[20] == 'Overlapping':
                        CSQ[0].extend(['uStop_gain to N-terminal extension'])
                        CSQ[1].extend(['N-terminal extension'])
                        Anno = uORF[1:3] + [uORF[4]] + uORF[17:19] + [uORF[19]] + uORF[20:-4] + uORF[-3:]
                        uORFAnnotations += [Anno]
                    elif codon > uSTOP and uORF[20] != 'Overlapping':
                        if codon == startPOS and uORF[20] != 'N-terminal extension':
                            CSQ[0].extend(['uStop_loss to N-terminal extension'])
                            CSQ[1].extend(['N-terminal extension'])
                            Anno = uORF[1:3] + [uORF[4]] + uORF[17:19] + [uORF[19] + " > " + NewUstopCodon] + uORF[20:-4] + uORF[-3:]
                            uORFAnnotations += [Anno]
                        elif codon > startPOS:
                            CSQ[0].extend(['uStop_loss to Overlapping'])
                            CSQ[1].extend(['decreased'])
                            Anno = uORF[1:3] + [uORF[4]] + uORF[17:19] + [uORF[19] + " > " + NewUstopCodon] + uORF[20:-4] + uORF[-3:]
                            uORFAnnotations += [Anno]
                        elif uORF[20] == 'Non-Overlapping':
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
            result.append(variant[-2].split('_') + variant[5:-4] + [variant[-4], variant[-1]] + [hit, CSQ[1][count]] + UTR[1:12] + UTR[14:] + uORFAnnotations[count])
            count += 1
    return result

def process_variants_spliceai_3(input_file_path, output_file_path, data_dir='~/.5ULTRA/data'):
    """Processes all variants and writes the results to the output file."""
    UTR_FILE_PATH = os.path.join(os.path.expanduser(data_dir), '5UTRs.tsv')
    UORF_FILE_PATH = os.path.join(os.path.expanduser(data_dir), 'uORFs.tsv')
    UTRs = load_tsv_data(UTR_FILE_PATH)
    uORFs = load_tsv_data(UORF_FILE_PATH)
    variants = load_tsv_data(input_file_path)
    utrs_by_transcript = defaultdict(list)
    uorfs_by_transcript = defaultdict(list)
    for UTR in UTRs[1:]:
        TRANSCRIPT = UTR[6]
        utrs_by_transcript[TRANSCRIPT].append(UTR)
    UTR_headers = UTRs[0]
    for uORF in uORFs[1:]:
        TRANSCRIPTS = uORF[5]
        uorfs_by_transcript[TRANSCRIPTS].append(uORF)
    uORF_headers = uORFs[0]
    with open(output_file_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        fields = variants[0][:-4] + [variants[0][-4], variants[0][-1]] + ['CSQ', 'translation'] + UTR_headers[1:12] + UTR_headers[14:] + uORF_headers[1:3] + [uORF_headers[4]] + uORF_headers[17:-4] + uORF_headers[-3:]
        writer.writerow(fields)
        for variant in variants[1:]:
            if ',' in variant[4]:
                continue
            processed_variant = process_variant_spliceai_3(variant, utrs_by_transcript, uorfs_by_transcript, data_dir)
            if processed_variant:
                writer.writerows(processed_variant)

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Detect variants.')
    parser.add_argument('input_file_path', type=str, help='Path to the filtered input file.')
    parser.add_argument('output_file_path', type=str, help='Path to the detection output file.')
    parser.add_argument('--data-dir', type=str, default='~/.5ULTRA/data', help='Path to the data directory.')
    args = parser.parse_args()

    process_variants_spliceai_3(args.input_file_path, args.output_file_path, data_dir=args.data_dir)

if __name__ == "__main__":
    main()
