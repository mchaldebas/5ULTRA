# 5ULTRA: A Pipeline for 5' UTR Variant Annotation and Scoring

**5ULTRA** is a computational pipeline designed to annotate and score genetic variants located in the 5′ untranslated regions (5′ UTRs) of genes. By focusing on upstream open reading frames (uORFs), Kozak sequences, and optional splicing sites (via SpliceAI), it provides detailed insights into how 5′ UTR variants can affect gene regulation, translation efficiency, and disease pathogenesis.

---

## Table of Contents
1. [Overview](#overview)  
2. [Features](#features)  
3. [Installation](#installation)  
4. [Usage](#usage)  
   - [Command-Line Options](#command-line-options)  
   - [Examples](#examples)  
5. [Pipeline Workflow](#pipeline-workflow)  
6. [Input and Output File Format](#input-and-output-file-format)  
7. [Reference](#reference)  
8. [Contributing](#contributing)  
9. [License](#license)  
10. [Contact](#contact)  

---

## Overview

The **5ULTRA** pipeline filters variants that reside within the 5′ UTR region and characterizes them for their potential impact on translation. Leveraging a machine learning model (random forest), the pipeline computes a comprehensive score, helping prioritize variants that could alter mRNA splicing or the regulation of the coding sequence (CDS).

**Key Highlights**:
- Pinpoints uORF changes (e.g., gain/loss of start or stop codons).
- Evaluates Kozak sequence disruptions.
- Integrates optional splicing analysis through [SpliceAI](https://github.com/Illumina/SpliceAI).
- Provides a single, interpretable score for rapid variant prioritization.

---

## Features

- **5′ UTR Filtering**: Retains only the variants that overlap with 5′ UTR regions.  
- **uORF & Kozak Annotation**: Identifies variants affecting uORFs (start/stop gain/loss) and Kozak sequences.  
- **Scoring**: Computes a comprehensive random forest–based score, summarizing the potential impact on translation.  
- **Optional Splicing Integration**: Includes [SpliceAI](https://github.com/Illumina/SpliceAI) predictions when `--splice` is enabled, offering insights into splice site disruptions within the 5′ UTR.

---

## Installation

### Prerequisites
- Python ≥ 3.6  
- [pip](https://pip.pypa.io/en/stable/)

### Steps

1.  **Clone the repository**  
    ```bash
    git clone https://github.com/mchaldebas/5ULTRA.git
    ```
2.  **Navigate into the 5ULTRA directory**  
    ```bash
    cd 5ULTRA
    ```
3.  **Install locally via pip**  
    ```bash
    pip install .
    ```
4.  **Test installation**
    ```bash
    5ULTRA -h
    ```
5.  **Download Required Data** (default path: ~/.5ULTRA/data)
    ```bash
    5ULTRA-download-data [--data-dir [path/to/data]]
    ```
    Specify --data-dir if you wish to store reference data in a non-default location.

**Note**: Ensure you have ~7 GB of available disk space to accommodate all reference data.
## Usage
Once installed, **5ULTRA** can be run directly as a command-line tool:
```bash
5ULTRA [-h] -I [input] [-O [output]] [--data-dir [path/to/data]] [--splice] [--full] [--mane]
```
### Options
- ```-h```: Show help message and exit
- ```-I```: Path to the input **VCF** or **TSV** file containing genetic variants.
- ```-O [output]```: Path for the output TSV file. Defaults to <input_file>.5ULTRA.tsv if not specified.
- ```--data-dir [path/to/data]```: Path to the data directory. Defaults to ~/.5ULTRA/data if not specified.
- ```--splice```: Enable SpliceAI processing and outputs on the variants with impact on splicing.
- ```--full```: Outputs a more detailed annotation (see Full Annotation columns under [Input and Output File Format](#input-and-output-file-format)).
-```--mane```: Outputs only the variants affecting MANE transcripts
### Examples
#### Basic Usage
```bash
5ULTRA tests/test-variants.tsv
```
This command reads test-variants.tsv, filters for 5′ UTR variants, annotates them, calculates scores, and writes the output to test-variants.5ULTRA.tsv.
#### Usage with all Arguments
```bash
5ULTRA  -I tests/test-variants.tsv \\
        -O tests/fully_annotated_variants.tsv \\
        --data-dir ~/.5ULTRA/data \\
        --splice \\
        --mane \\
        --full
```
This command uses custom data paths, analyse only splicing variants of MANE transcripts, and produces a fully annotated output with additional columns.
## Input and Output File Format

- ***Input***: **VCF** or **TSV** file with genetic variants. (#CHROM, POS, ID, REF, ALT)
- ***Output***: TSV file containing annotated and scored variants.
    - **CHROM, POS, ID, REF, ALT**: Same as input.
    - **CSQ**: Type of variant
    - **Translation**: Increased, Decreased, or N-terminal Extension
    - **5ULTRA_Score**: Prioritization metric
    - **SpliceAI** (optional): SpliceAI scores for splicing variants that may affect 5’UTR sequence
    - **Splicing_CSQ** (optional): Type of splicing variant
    - **GENE**: Gene Symbol
    - **TRANSCRIPT**: Ensembl transcript ID 

- ***Full Annotation*** 
When --full is specified, additional columns are appended:
    - **MANE**: NCBI transcript ID if applicable (e.g., NM_123456789.1)
    - **5UTR_START / 5UTR_END**: Genomic coordinates of the 5′ UTR.
    - **STRAND**: DNA strand (+ or -).
    - **5UTR_LENGTH**: Length of the 5’UTR.
    - **START_EXON**: The exon number where the CDS starts.
    - **mKOZAK / mKOZAK_STRENGTH**: Kozak sequence/context strength around the CDS start.
    - **uORF_count**: Total number of uORFs in the transcript.
    - **Overlapping_count, Nterminal_count, NonOverlapping_count**: Counts of specific uORF types.
    - **uORF_START / uORF_END**: Genomic coordinates of the uORF.
    - **Ribo_seq**: Indicates evidence of translation  (*True*, *False* or *New uORF*)
    - **uSTART_mSTART_DIST / uSTART_CAP_DIST**: Distance from CDS start and 5′ cap, respectively.
    - **uSTOP_CODON**: Type of stop codon (TAA, TGA, or TAG).
    - **uORF_TYPE**: uORF type (Non-overlapping, Overlapping, N-terminal extension).
    - **uKOZAK / uKOZAK_STRENGTH**: Kozak sequence/context strength around the uORF start.
    - **uORF_LENGTH / uORF_AA_LENGTH**: Length of the uORF in nucleotides and amino acids.
    - **uORF_rank**: Relative rank based on proximity to the CDS.
    - **uSTART_PHYLOP / uSTART_PHASTCONS**: Conservation scores (PhyloP, PhastCons).
    - **pLI / LOEUF**: Gene-level intolerance metrics.


## Reference

***Chaldebas M. et al.*** (2025). *Genome-wide detection of human 5’UTR candidate variants.*

## Contributing

Contributions are welcome! Please submit pull requests or open issues on the GitHub repository.

## License

This project is licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International ([CC BY-NC-ND 4.0](LICENSE)). See the LICENSE file for details.

## Contact
> **Developer:** Matthieu Chaldebas, Ph.D. candidate

> **Email:** mchaldebas@rockefeller.edu

> **Laboratory:** St. Giles Laboratory of Human Genetics of Infectious Diseases

> **Institution:** The Rockefeller University, New York, NY, USA

For any questions, feel free to email directly.
