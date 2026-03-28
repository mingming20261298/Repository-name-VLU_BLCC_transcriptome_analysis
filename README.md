# VLU BLCC Transcriptome Analysis

## Overview
This repository contains analysis scripts and results for a secondary transcriptomic study of venous leg ulcers (VLUs) using dataset GSE84571. The study applies a paired difference-in-differences (DiD) framework to investigate early molecular responses associated with BLCC (Apligraf) treatment.

## Study Design
- Dataset: GEO accession GSE84571  
- Platform: Affymetrix GPL570  
- Design: Paired samples (week 0 vs. week 1)  
- Groups:
  - Compression therapy (control)
  - BLCC + compression (treatment)

## Analysis Workflow
The analysis includes the following steps:

1. Preprocessing
   - RMA normalization of CEL files
   - Probe-to-gene mapping (GPL570 annotation)

2. Differential Expression
   - limma with paired design (patient blocking)
   - Difference-in-Differences (DiD) model

3. Pathway Analysis
   - GSEA (fgsea)
   - GSVA pathway scoring

4. Microenvironment Estimation
   - xCell / MCP-counter

5. Signature Construction
   - Weighted gene score
   - Penalized regression (LASSO / Elastic Net)

## Repository Structure
.
├── scripts/
├── data/
├── results/
├── figures/

## Data Availability
The dataset is publicly available from GEO:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84571

## Reproducibility
All analyses were performed in R using packages including limma, GSVA, fgsea, affy/oligo, xCell, and MCP-counter.

## Code Availability
Analysis scripts and processed data supporting the findings are available in this repository.

Archived version:
https://doi.org/10.5281/zenodo.19042183

## Notes
This repository provides the main analytical steps and results used in the manuscript. Minor differences may occur due to environment variations.
