# Co-expression Network Analysis of APOE Genotypes in Alzheimer's Disease



## Project Overview

This project investigates the impact of APOE genotypes on gene co-expression networks in Alzheimer's disease (AD) using RNA-sequencing data from two major cohorts: the Religious Orders Study and Memory and Aging Project (ROSMAP) and the Mount Sinai Brain Bank (MSBB). The primary goal is to identify co-expression patterns that are altered across different APOE genotypes (E2E3, E3E3, and E3E4), providing insights into genotype-specific molecular mechanisms underlying AD pathology. The analysis workflow is designed to be fully reproducible, from raw data processing to the identification of genotype-specific networks and their functional enrichment.



## Key Features

- **Cross-Cohort Analysis:** Integration of ROSMAP and MSBB datasets to identify robust and reproducible co-expression patterns.

- **Genotype-Specific Networks:** Construction of co-expression networks independently for APOE-E2E3, E3E3, and E3E4 carriers.

- **Pattern Identification:** Development of a robust methodology to filter for monotonic co-expression patterns, such as a decreasing or increasing correlation strength across the E2E3 -> E3E3 -> E3E4 axis.

- **Network Exploration:** Detailed analysis of network connections for key genes, including APOE and other AD-related genes.

- **Functional Enrichment:** Gene Ontology (GO) and KEGG pathway enrichment analysis to determine the biological functions associated with the identified patterns.



## Repository Structure

The project is organized into a modular structure to ensure clarity and reproducibility:



- `data/`: Contains all raw and reference data files. The `raw/` subdirectory (with large files like MSBB raw counts) is excluded from version control via `.gitignore`.

- `processed\_data/`: Stores the cleaned and normalized R objects (`.rds` files) ready for analysis.

- `output/`: The destination for all final analysis results, including Excel reports and plots.

- `patterns/`: Holds the data frames containing the identified co-expression patterns.

- `scripts/`: Houses all R scripts, organized by function.

&nbsp;   - `\*\_data\_processing.R`: Scripts for handling raw data, quality control, and normalization.

&nbsp;   - `\*\_functions.R`: A library of helper functions for correlation, network, and enrichment analysis.

&nbsp;   - `main\_analysis.R`: The central script that orchestrates the entire workflow.



## Prerequisites

To run this analysis, you need to have R installed along with the following packages. You can install them using the commands below:

```R

\# Install CRAN packages

install.packages(c("dplyr", "tidyverse", "data.table", "ggplot2", "writexl", "Hmisc", "foreach", "doParallel"))



\# Install Bioconductor packages

if (!requireNamespace("BiocManager", quietly = TRUE))

&nbsp;   install.packages("BiocManager")

BiocManager::install(c("DESeq2", "limma", "edgeR", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "biomaRt"))

