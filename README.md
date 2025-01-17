# 5ULTRA: A Pipeline for 5' UTR Variant Annotation and Scoring

## Overview

5ULTRA is a computational pipeline designed to annotate and score genetic variants located in the 5' untranslated regions (5' UTRs) of genes. The pipeline focuses on identifying variants that may affect upstream open reading frames (uORFs), Kozak sequences, and splicing sites, which can have significant implications in gene regulation and disease.

## Features

- 5' UTR Filtering: Filters input variant files to retain only those variants located within 5' UTR regions.
- Variant Annotation: Detects and Annotates variants that impact uORFs, including uStart gains/losses, uStop gains/losses, and Kozak sequence alterations.
- Scoring: Calculates a comprehensive score for each variant based on a random forest model, aiding in the interpretation of potential effect on CDS translation.
- SpliceAI Integration: Optionally incorporates SpliceAI predictions to assess the impact of variants on splicing.

## Pipeline Workflow

1. Input Preparation: Accepts VCF or TSV files containing genetic variants.
2. 5' UTR Filtering: Filter-input.py filters for variants within 5' UTR regions based on intervals.
3. Variant Detection:
    - Detection.py annotates variants affecting uORFs and Kozak sequences.
    - If SpliceAI is enabled, executes Spliceai-Main.sh for splicing impact analysis.
4. Scoring: Applies Score.py to compute a variant score using a pre-trained machine learning model.
5. Output Generation: Produces a TSV file with annotated and scored variants.

## Installation

### Prerequisites
- Python >=3.6
- pip
### Steps
1. Pip installation:
```
pip install fiveULTRA
```
2. Download data (default path: ~/.5ULTRA/data):
```
5ULTRA-download-data [--data-dir [path/to/data]]
```
## Usage
```
5ULTRA [-h] -I [input] [-O [output]] [--data-dir [path/to/data]] [--splice] [--full]
```
### Options
- ```-h```: Show help message and exit
- ```-I```: Path to the input VCF or TSV file containing genetic variants.
- ```-O [output]```: Path for the output TSV file. Defaults to <input_file>.5ULTRA.tsv if not specified.
- ```--data-dir [path/to/data]```: 
- ```--splice```: Enable SpliceAI processing for splicing impact analysis.
- ```--full```: Get the Full annotation in the output.
### Examples
#### Basic Usage
```
python 5ULTRA.py test-variants.tsv
```
#### Usage with all Arguments
```
python 5ULTRA.py -I test-variants.tsv -O fully_annotated_variants.tsv --data-dir ~/.5ULTRA/data --splice --full 
```
## Input and Output File Format

- ***Input***: VCF or TSV file with genetic variants.
- ***Output***: TSV file containing annotated and scored variants.
    - CHROM, POS, ID, REF, ALT (same as input)
    - CSQ: Type of variant
    - Translation: Increased, Decreased, or N-terminal Extension
    - 5ULTRA_Score: Prioritization metric
    - [SpliceAI]: SpliceAI scores for splicing variants that may affect 5’UTR sequence
    - [Splicing_CSQ]: Type of splicing variant
    - GENE: Gene Symbol
    - TRANSCRIPT: Ensembl transcript ID 

- ***Full Annotation***: Additional columns in the output when --full specified
    - MANE: NCBI transcript ID if applicable (e.g., NM_123456789.1)
    - 5UTR_START: Genomic position of the 5’UTR start.
    - 5UTR_END: Genomic position of the 5’UTR end.
    - STRAND: DNA strand (+ or -).
    - 5UTR_LENGTH: Length of the 5’UTR.
    - START_EXON: CDS start exon position.
    - mKOZAK: Nucleotide sequence -4 to +5 around the CDS start.
    - mKOZAK_STRENGTH: CDS Kozak context strength (Weak, Adequate, Strong, or NA).
    - uORF_count: Total number of uORFs in the transcript.
    - Overlapping_count: Number of overlapping uORFs.
    - Nterminal_count: Number of N-terminal extension uORFs.
    - NonOverlapping_count: Number of non-overlapping uORFs.
    - uORF_START: Genomic position of the uORF start.
    - uORF_END: Genomic position of the uORF end.
    - Ribo_seq: Evidence of translation (True, False or New uORF)
    - uSTART_mSTART_DIST: uORF start distance from the CDS start.
    - uSTART_CAP_DIST: uORF start distance from the 5’UTR cap.
    - uSTOP_CODON: Type of stop codon (TAA, TGA, or TAG).
    - uORF_TYPE: uORF type (Non-overlapping, Overlapping, N-terminal extension).
    - uKOZAK: Nucleotide sequence -4 to +5 around the uORF start.
    - uKOZAK_STRENGTH: uORF Kozak context strength (Weak, Adequate, Strong, or NA).
    - uORF_LENGTH: Length of the uORF.
    - uORF_AA_LENGTH: Amino Acid length of the uORF.
    - uORF_rank: Rank of the uORF based on proximity to the CDS start.
    - uSTART_PHYLOP: Mean conservation score of the uORF start (PhyloP).
    - uSTART_PHASTCONS: Mean conservation score of the uORF start (PhastCons).
    - pLI: Gene probability of being loss-of-function intolerant.
    - LOEUF: Gene loss-of-function observed/expected upper bound fraction.


## Reference

***Chaldebas M. et al.*** Genome-wide detection of human 5’UTR candidate variants. 2025.

## Contributing

Contributions are welcome! Please submit pull requests or open issues on the GitHub repository.

## License

This project is licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International ([CC BY-NC-ND 4.0](LICENSE)). See the LICENSE file for details.

## Contact
> **Developer:** Matthieu Chaldebas, Ph.D. candidate

> **Email:** mchaldebas@rockefeller.edu

> **Laboratory:** St. Giles Laboratory of Human Genetics of Infectious Diseases

> **Institution:** The Rockefeller University, New York, NY, USA
