# File: msbb_data_processing.R
# This script processes and cleans the MSBB transcriptome data and its associated metadata.
# It handles data import from raw featureCounts files, metadata integration,
# filtering, normalization, and quality control (PCA) to prepare the final
# data for co-expression network analysis.

# 1. Load Libraries
library(dplyr)
library(tidyverse)
library(tximport)
library(readr)
library(DESeq2)
library(limma)
library(edgeR)
library(ggplot2) # For PCA plotting

# 2. Set Working Directory
# Please set this to your main project directory.
setwd('path/APOE_AD_Stud/')

# 3. Load Raw Count Data and Prepare Matrix
# Set the directory containing the featureCounts files.
dir <- "data/raw/MSBB/"
files <- list.files(path = dir, pattern = "\\.txt$", full.names = TRUE)

# Read the gene IDs from the first file.
sample1 <- read_tsv(files[1], comment = "#")
gene_ids <- sample1$Geneid

# Read counts from all files, assuming the 7th column contains the count data.
read_counts <- function(file) {
  df <- read_tsv(file, comment = "#", col_names = FALSE)
  df %>% select(X7)
}
counts_list <- lapply(files, read_counts)
counts_matrix <- do.call(cbind, counts_list)
counts_matrix <- counts_matrix[-1, ]

# Set gene IDs as row names and sample names from file names.
rownames(counts_matrix) <- gene_ids
colnames(counts_matrix) <- basename(files)
colnames(counts_matrix) <- gsub("_featurecounts\\.txt", "", colnames(counts_matrix))

# Convert the matrix to numeric values.
counts_matrix <- as.data.frame(lapply(counts_matrix, as.numeric))
rownames(counts_matrix) <- gene_ids
counts_matrix <- as.matrix(counts_matrix)

# 4. Prepare Metadata
# Load and filter metadata to match samples in the count matrix.
biospecimen_metadata <- read.csv('data/reference/FP_metadata.csv')
biospecimen_metadata <- biospecimen_metadata %>% filter(specimenID %in% colnames(counts_matrix))

# 5. Filter for Specific APOE Genotypes and CERAD Score
# Filter for CERAD score 1 and remove the E2E2 genotype.
meta_data_final <- biospecimen_metadata %>% filter(CERAD == 1 & apoeGenotype != 22)
meta_data_final <- na.omit(meta_data_final)

# Update APOE genotype factor levels for consistency.
apoe_categ <- c("23" = "E2E3", "33" = "E3E3", "34" = "E3E4")
meta_data_final <- meta_data_final %>% mutate(
  apoeGenotype = factor(apoeGenotype, levels = names(apoe_categ), labels = apoe_categ)
)
names(meta_data_final)[names(meta_data_final) == 'apoeGenotype'] <- 'apoe_genotype'

# Clean gene IDs in the count matrix to remove version numbers.
gene_ids_cleaned <- sub("\\..*$", "", rownames(counts_matrix))
counts_matrix_cleaned <- counts_matrix
rownames(counts_matrix_cleaned) <- gene_ids_cleaned

# 6. Filter Low-Expressed Genes and Align Data
# Filter count data to match the final metadata.
counts_final <- counts_matrix_cleaned[, colnames(counts_matrix_cleaned) %in% meta_data_final$specimenID]
counts_final <- as.matrix(counts_final)

# Filter low-expressed genes using edgeR.
dge <- DGEList(counts = counts_final)
keep <- filterByExpr(dge)
dge <- dge[keep, keep.lib.size = FALSE]
counts_filtered <- counts_final[keep, ]
counts_filtered <- as.matrix(counts_filtered)

# 7. Quality Control with PCA
# Create a DESeq2 object for VST normalization and PCA.
meta_data_final$sex <- factor(meta_data_final$sex)
meta_data_final$ageDeath <- as.numeric(gsub('90\\+','90', meta_data_final$ageDeath))

dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = meta_data_final, design = ~1)
dds <- DESeq(dds)
vst <- vst(dds)

# Perform PCA and visualize clustering by APOE genotype.
pca.norm <- prcomp(t(assay(vst)), scale = TRUE)
pca.var.norm <- pca.norm$sdev^2
pca.var.per.norm <- round(pca.var.norm/sum(pca.var.norm)*100,1)
pca.data.norm <- data.frame(Sample=rownames(pca.norm$x),
                            X=pca.norm$x[,1],
                            Y=pca.norm$x[,2],
                            status = meta_data_final$apoe_genotype)

ggplot(data=pca.data.norm, aes(x=X, y=Y, colour=factor(status))) +
  geom_point() +
  labs(x = paste('PCA1-', pca.var.per.norm[1], '%', sep = ''),
       y = paste('PCA2-', pca.var.per.norm[2], '%', sep = ''),
       colour = 'status') +
  theme_minimal()

# 8. Save Final Processed Data
# Save the final filtered data and metadata in an RDS file.
msbb_processed_data <- list(counts = counts_filtered, meta_data = meta_data_final)
saveRDS(msbb_processed_data, 'output/msbb_processed_data.rds')

