# Genome Analysis: Stop Codons

Scripts to process mitochondrial DNA (mtDNA) stop codon data across various taxa. This code was used for analysis in the review article on translation termination in human mitochondria.

## Requirements

- Python > 3.8 with Jupyter notebook support
- MacOS/Linux system

## Repository Structure

- `src/`: Contains utilities for GenBank file processing.
- `processing/`: Contains the Jupyter notebook (`stop_codon_analysis.ipynb`) with analytics on stop codons.

## Getting Started

To reproduce the analysis, follow these steps:

### Clone the repository

```bash
git clone git@github.com:dm-shr/genome_analysis_stop_codons.git
cd genome_analysis_stop_codons
```

### Run the Makefile

This will download the annotated RefSeq mitochondrial genomes in GenBank format and process them to prepare for the analysis.

```bash
make all
```

### Analysis Notebook

Next, follow the cells in the processing/stop_codon_analysis.ipynb notebook to reproduce the analysis.

### Cleaning Up
To delete the data processed during the analysis, you can use the make clean command.

```bash
make clean
```
