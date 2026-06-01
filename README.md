# Differential expression analysis associated with lactic acid in murine macrophages during M1 activation

This repository contains the scripts, metadata, summarized results, figures, and reproducibility files for the final Transcriptomics project.

## Overview

This project implements a reproducible paired-end RNA-seq analysis workflow using public data from the GSE115354 series. The workflow includes metadata filtering, FASTQ download from ENA, quality control, read cleaning, alignment with STAR, gene-level quantification with featureCounts, differential expression analysis with DESeq2, and functional enrichment analysis with GO Biological Process.

## Repository structure

- `01_metadata/`: processed metadata and SRR accession lists.
- `03_scripts/`: scripts used throughout the analysis workflow.
- `04_results/`: summarized tables and main analysis results.
- `05_figures/`: figures generated for the report.
- `reproducibility/`: software versions, R package versions, and computational environment documentation.
- `Reporte_proyecto.qmd`: source file for the final report.
- `Reporte_proyecto.pdf`: final report in PDF format.

## Data not included

Large files such as FASTQ files, BAM files, reference genome files, annotation files, and STAR indexes are not included in this repository due to size limitations. These files can be regenerated using the SRR accessions and the scripts provided in the repository.

## Dataset and reference files

- Dataset: GSE115354
- FASTQ source: ENA
- Organism: *Mus musculus*
- Reference genome: mm39
- Annotation: GENCODE vM36

## Main tools

The main tools used in this workflow were:

- FastQC and MultiQC for quality control.
- fastp for read cleaning.
- STAR for genome alignment.
- featureCounts for gene-level quantification.
- DESeq2 for differential expression analysis.
- GO Biological Process for functional enrichment analysis.

## Reproducibility

Software versions, R package versions, and computational environment information are documented in the `reproducibility/` directory. The scripts are numbered according to the order of execution in the workflow.
