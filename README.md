# 5ULTRA: A Pipeline for 5' UTR Variant Annotation and Scoring

## Overview

5ULTRA is a computational pipeline designed to annotate and score genetic variants located in the 5' untranslated regions (5' UTRs) of genes. The pipeline focuses on identifying variants that may affect upstream open reading frames (uORFs), Kozak sequences, and splicing sites, which can have significant implications in gene regulation and disease.

## Features

- 5' UTR Filtering: Filters input variant files to retain only those variants located within 5' UTR regions.
- Variant Annotation: Detects and annotates variants that impact uORFs, including uStart gains/losses, uStop gains/losses, and Kozak sequence alterations.
- SpliceAI Integration: Optionally incorporates SpliceAI predictions to assess the impact of variants on splicing.
- Scoring: Calculates a comprehensive score for each variant based on a reandom forest model, aiding in the interpretation of potential effect on CDS translation.

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
- Operating System: Linux or macOS
- Software:
    - Python 3.x
- Python Packages:
    - pysam
    - pandas
    - numpy
    - joblib
    - sklearn
### Steps
1. Clone the Repository:
```
git clone https://github.com/mchaldebas/5ULTRA.git
cd 5ULTRA
```
2. Install Python Dependencies:
```
pip install pysam pandas numpy joblib sklearn
```
3. Ensure Data Files Are in Place:
    - Place 5UTRs.intervals.bed, 5UTRs.tsv, and uORFs.tsv in the ./data directory.
    - Ensure PhyloP and PhastCons score files are available and indexed with tabix.

## Usage
```
python 5ULTRA.py [--splice] <input_file> [output_file]
```
### Options
- --splice: Enable SpliceAI processing for splicing impact analysis.
- <input_file>: Path to the input VCF or TSV file containing genetic variants.
- [output_file]: (Optional) Path for the output TSV file. Defaults to <input_file>.5ULTRA.tsv if not specified.
### Examples
#### Basic Usage
```
python 5ULTRA.py input_variants.vcf
```
#### Usage with SpliceAI Integration to look only at splicing variants
```
python 5ULTRA.py --splice input_variants.vcf annotated_variants.tsv
```
## Input and Output File Format

- ***Input***: VCF or TSV file with genetic variants.
- ***Output***: TSV file containing annotated and scored variants.

    - CHROM, POS, ID, REF, ALT (exactly the same as input)
    - CSQ
    - GENE
    - CDS translation
    - 5ULTRA_Score
- ***Full Annotation***:
    - 
    - 
    - 
    - 
    - 
    - 
    - 
    - 
    - 

## Reference

***Chaldebas M. et al.*** Genome-wide detection of human 5’UTR candidate variants. 2025.

## Contributing

Contributions are welcome! Please submit pull requests or open issues on the GitHub repository.

## License

This project is licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.

## Contact
> **Developer:** Matthieu Chaldebas, Ph.D. candidate

> **Email:** mchaldebas@rockefeller.edu

> **Laboratory:** St. Giles Laboratory of Human Genetics of Infectious Diseases

> **Institution:** The Rockefeller University, New York, NY, USA

## Disclaimer
This tool is intended for research purposes only and is not certified for clinical use.