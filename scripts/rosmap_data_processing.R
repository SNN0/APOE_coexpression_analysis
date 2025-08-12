# File: rosmap_data_processing.R
# This script processes and cleans the ROSMAP transcriptome data and its associated metadata.
# It performs data loading, merging, filtering, quality control (PCA), and normalization
# to prepare the final data for co-expression network analysis.

# 1. Load Libraries
library(dplyr)
library(tidyverse)
library(DESeq2)
library(limma)
library(edgeR)
library(ggplot2) # For PCA plotting

# 2. Set Working Directory
# Please set this to your main project directory.
setwd('path/APOE_AD_Stud/')

# 3. Load Data
# Assuming 'rosmapallBusra.RData' contains the initial 'rosmapcountdata' object
# and the other files are in the same directory.
load('data/raw/ROSMAP/rosmapallBusra.RData')

# Load raw count data and clinical metadata.
count_data_main <- read.csv('data/raw/ROSMAP/count_data.csv', stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE, header = TRUE)
biospecimen_metadata <- read.csv('data/raw/ROSMAP/ROSMAP_biospecimen_metadata.csv')
main_meta_data <- read.csv('data/raw/ROSMAP/ROSMAP_clinical.csv')

# 4. Prepare Count Data
# Clean and align column names of the count data.
colnames(count_data_main)[1] <- 'bos'
count_data_main <- count_data_main %>% select(-bos)
col.names.count <- colnames(rosmapcountdata)
colnames(count_data_main) <- col.names.count

# Extract gene symbols, assuming first two columns contain gene info.
gene_symbol <- count_data_main[, c(1, 2)]
# Filter out the gene info columns to get the count matrix and remove an outlier sample.
counts <- count_data_main[, 5:641]
counts <- counts %>% select(-'764_130520')

# 5. Prepare Metadata
# Filter biospecimen metadata to match samples in count data.
speciman_ids <- colnames(counts)
distinct_biospecimen_data <- biospecimen_metadata %>% distinct(specimenID, .keep_all = TRUE)
biospecimen_data_sub <- distinct_biospecimen_data[distinct_biospecimen_data$specimenID %in% speciman_ids, ]

# Filter clinical metadata to match individuals in biospecimen data.
info_subset <- c('projid', 'msex', 'apoe_genotype', 'age_death', 'pmi', 'ceradsc', 'individualID')
main_meta_data_subset <- main_meta_data[, info_subset]
individualIDs <- biospecimen_data_sub$individualID
meta_data <- main_meta_data_subset[main_meta_data_subset$individualID %in% individualIDs, ]

# Merge clinical and biospecimen data.
meta_data <- merge(meta_data, biospecimen_data_sub[,c(1,2)], by = 'individualID', all.x = TRUE)

# 6. Metadata Annotation and Cleaning
# Annotate AD status based on ceradsc scores.
meta_data <- meta_data %>% mutate(status = case_when(
  ceradsc %in% c(1, 2) ~ 'AD',
  ceradsc %in% c(3, 4) ~ 'Control',
  TRUE ~ NA_character_
))

# Define APOE and sex categories and convert to factors.
apoe_categ <- c("23" = "E2E3", "33" = "E3E3", "34" = "E3E4", "22" = "E2E2", "24" = "E2E4", "44" = "E4E4")
msex_categ <- c('0' = 'Female', '1' = 'Male')
meta_data <- meta_data %>% mutate(
  apoe_genotype = factor(apoe_genotype, levels = names(apoe_categ), labels = apoe_categ),
  msex = factor(msex, levels = names(msex_categ), labels = msex_categ)
)
meta_data$age_death <- as.numeric(gsub('90\\+','90', meta_data$age_death))

# Re-order metadata to match the count matrix and omit rows with NA values.
meta_data <- meta_data[order(meta_data$specimenID), ]
meta_data <- na.omit(meta_data)

# 7. Filter Low-Expressed Genes and Align Data
# Align counts and metadata, then filter low-expressed genes.
counts <- counts %>% select(meta_data$specimenID)
counts <- as.matrix(counts)

dge <- DGEList(counts = counts, group = meta_data$status)
keep <- filterByExpr(dge, group = meta_data$status)
counts_filtered <- counts[keep, ]
counts_filtered <- as.matrix(counts_filtered)

# 8. Quality Control with PCA and Outlier Removal
# First, remove gene info from the count matrix rows for DESeq2 compatibility.
gene_info <- count_data_main[, 1:2]
gene_info_filtered <- gene_info[keep,]
rownames(counts_filtered) <- gene_info_filtered$ensembl_gene_id

# Create a DESeq2 object for VST normalization and PCA.
dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = meta_data, design = ~status + msex + age_death + pmi)
dds <- DESeq(dds)
vst <- vst(dds)

# Perform PCA and identify outlier `380_120503`.
pca.norm <- prcomp(t(assay(vst)), scale = TRUE)
pca.data.norm <- data.frame(Sample=rownames(pca.norm$x),
                            X=pca.norm$x[,1],
                            Y=pca.norm$x[,2],
                            status =meta_data$msex)
# Visualization of the PCA plot would show the outlier.
ggplot(data=pca.data.norm, aes(x=X, y=Y, colour=factor(status))) +
  geom_point() +
  labs(x = 'PCA1', y = 'PCA2', colour = 'msex') +
  theme_minimal()

# Remove the identified outlier sample.
meta_data <- meta_data %>% filter(specimenID != '380_120503')
counts_filtered <- counts_filtered %>% as.data.frame() %>% select(meta_data$specimenID)
counts_filtered <- as.matrix(counts_filtered)

# 9. Final Filtering for Co-expression Analysis
# Keep only samples with CERAD score 1 and specific APOE genotypes.
meta_data_final <- meta_data %>% filter(ceradsc == 1) %>% 
  filter(apoe_genotype %in% c('E2E3', 'E3E3', 'E3E4'))
meta_data_final$apoe_genotype <- factor(meta_data_final$apoe_genotype, levels = c('E2E3', 'E3E3', 'E3E4'))

# Filter count data to match the final metadata.
counts_final <- counts_filtered %>% as.data.frame() %>% select(meta_data_final$specimenID)
counts_final <- as.matrix(counts_final)

# 10. Save Final Processed Data
# Save the final filtered data and metadata in an RDS file.
rosmap_processed_data <- list(counts = counts_final, meta_data = meta_data_final, gene_symbols = gene_info_filtered)
saveRDS(rosmap_processed_data, 'output/rosmap_processed_data.rds')
