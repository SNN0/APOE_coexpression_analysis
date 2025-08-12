
# This script orchestrates the entire co-expression network analysis workflow,
# from loading the processed data to performing network and enrichment analyses.

# 1. Setup Environment
# Set the main project directory.
setwd('path/APOE_AD_Stud/')

# Load all required libraries.
library(dplyr)
library(tidyverse)
library(DESeq2)
library(data.table)
library(biomaRt)
library(writexl)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)

# Source all helper functions from their respective files.
# Make sure these files are in the same directory.
source('scripts/rosmap_data_processing.R')
source('scripts/msbb_data_processing.R')
source('scripts/correlation_functions.R')
source('scripts/network_functions.R')
source('scripts/enrichment_functions.R')

# 2. Load Processed Data
# Make sure the data processing scripts have been run previously to generate these files.
rosmap_data <- readRDS('output/rosmap_processed_data.rds')
msbb_data <- readRDS('output/msbb_processed_data.rds')

# Extract data and metadata for convenience.
rosmap_counts <- rosmap_data$counts
rosmap_meta <- rosmap_data$meta_data
msbb_counts <- msbb_data$counts
msbb_meta <- msbb_data$meta_data
msbb_meta$apoe_genotype <- factor(msbb_meta$apoe_genotype)

# 3. Perform Normalization (VST) for Co-expression Analysis
# ROSMAP data VST transformation
dds_rosmap <- DESeqDataSetFromMatrix(countData = rosmap_counts, colData = rosmap_meta, design = ~1)
dds_rosmap <- DESeq(dds_rosmap)
norm_counts_rosmap <- assay(vst(dds_rosmap)) %>% t()

# MSBB data VST transformation
dds_msbb <- DESeqDataSetFromMatrix(countData = msbb_counts, colData = msbb_meta, design = ~1)
dds_msbb <- DESeq(dds_msbb)
norm_counts_msbb <- assay(vst(dds_msbb)) %>% t()

# 4. Correlation Analysis and Pattern Identification
genotypes <- c('E2E3', 'E3E3', 'E3E4')

# ROSMAP: Calculate correlations and filter patterns
correlations_rosmap <- calculate_correlations(norm_counts_rosmap, rosmap_meta, genotypes)
signi_pairs_rosmap_e2e3 <- filter_significant_pairs(correlations_rosmap, 'E2E3', threshold_corr = 0.95, threshold_p = 0.05)
pattern_rosmap_dec <- filter_correlation_patterns_E2E3_ref(signi_pairs_rosmap_e2e3, correlations_rosmap, genotypes, threshold_diff = 0.25)
signi_pairs_rosmap_e3e4 <- filter_significant_pairs(correlations_rosmap, 'E3E4', threshold_corr = 0.80, threshold_p = 0.05)
pattern_rosmap_inc <- filter_correlation_patterns_E3E4_ref(signi_pairs_rosmap_e3e4, correlations_rosmap, genotypes, threshold_diff = 0.25)

# Save ROSMAP patterns
saveRDS(pattern_rosmap_dec, 'pattern/Rosmap_0.95cor_decreasing.rds')
saveRDS(pattern_rosmap_inc, 'pattern/Rosmap_0.80cor_increasing.rds')

# MSBB: Calculate correlations and filter patterns
correlations_msbb <- calculate_correlations(norm_counts_msbb, msbb_meta, genotypes)
signi_pairs_msbb_e2e3 <- filter_significant_pairs(correlations_msbb, 'E2E3', threshold_corr = 0.70, threshold_p = 1)
pattern_msbb_dec <- filter_correlation_patterns_E2E3_ref(signi_pairs_msbb_e2e3, correlations_msbb, genotypes, threshold_diff = 0.1)
signi_pairs_msbb_e3e4 <- filter_significant_pairs(correlations_msbb, 'E3E4', threshold_corr = 0.70, threshold_p = 1)
pattern_msbb_inc <- filter_correlation_patterns_E3E4_ref(signi_pairs_msbb_e3e4, correlations_msbb, genotypes, threshold_diff = 0.1)

# Save MSBB patterns
saveRDS(pattern_msbb_dec, 'pattern/msbb_0.95cor_decreasing.rds')


# 5. Cross-Dataset Comparison
# Combine positive and negative decreasing patterns for intersection analysis.
rosmap_trans_dec <- rbind(pattern_rosmap_dec$positive_cor_decreasing, pattern_rosmap_dec$negative_cor_decreasing)
msbb_trans_dec <- rbind(pattern_msbb_dec$positive_cor_decreasing, pattern_msbb_dec$negative_cor_decreasing)
rosmap_trans_inc <- rbind(pattern_rosmap_inc$positive_cor_increasing, pattern_rosmap_inc$negative_cor_increasing)
msbb_trans_inc <- rbind(pattern_msbb_inc$positive_cor_increasing, pattern_msbb_inc$negative_cor_increasing)

# Find common unique genes (network-based)
rosmap_trans_dec_network <- read.csv('data/reference/rosmap_095_Decreasing_Network.csv')
msbb_trans_dec_network <- read.csv('None')
msbb_trans_dec_network$gene1 <- sub("\\..*", "", msbb_trans_dec_network$gene1)
msbb_trans_dec_network$gene2 <- sub("\\..*", "", msbb_trans_dec_network$gene2)
rosmap_trans_dec_network_unique <- unique(c(rosmap_trans_dec_network$gene1, rosmap_trans_dec_network$gene2))
msbb_trans_dec_network_unique <- unique(c(msbb_trans_dec_network$gene1, msbb_trans_dec_network$gene2))
common_genes_dec <- intersect(rosmap_trans_dec_network_unique, msbb_trans_dec_network_unique)

# Find common pairs (non-network based)
common_pairs_dec <- find_common_pairs(rosmap_trans_dec, msbb_trans_dec)
common_pairs_inc <- find_common_pairs(rosmap_trans_inc, msbb_trans_inc)

# 6. Enrichment Analysis
# Perform GO and KEGG enrichment on the common gene list.
go_bp_results <- enrichment_go(common_genes_dec, 'BP', 'Common Decreasing Genes GO-BP')
go_cc_results <- enrichment_go(common_genes_dec, 'CC', 'Common Decreasing Genes GO-CC')
go_mf_results <- enrichment_go(common_genes_dec, 'MF', 'Common Decreasing Genes GO-MF')
kegg_results <- enrichment_kegg(common_genes_dec, 'Common Decreasing Genes KEGG')

# Combine and visualize results .
all_enrichment_results <- rbind(
  go_bp_results@result %>% mutate(Category = "BP"),
  go_cc_results@result %>% mutate(Category = "CC"),
  go_mf_results@result %>% mutate(Category = "MF"),
  kegg_results@result %>% mutate(Category = "KEGG")
)
write_xlsx(all_enrichment_results, 'output/common_decreasing_enrichment.xlsx')

# 7. Network Exploration 
# APOE-related connections
rosmap_trans_all_pattern <- rbind(rosmap_trans_dec, rosmap_trans_inc)
apoE_connections <- find_all_connections(rosmap_trans_all_pattern, 'ENSG00000130203', degree_limit = 2)
write_xlsx(apoE_connections, 'output/APOE_connections_rosmap.xlsx')
# Plotting specific gene-pair correlations (example)
# plot_gene_correlation(rosmap_meta, norm_counts_rosmap, 'ENSG00000130203', 'ENSG00000087258')

# Connections for a list of AD-related genes
AD_related_gene_list <- read.csv('AD_related_gene_list.yeni.csv')
rosmap_dec_adrelated <- intersect(rosmap_trans_dec_unique, AD_related_gene_list$ENSEMBL)
ad_gene_network <- find_all_connections_combined(rosmap_trans_dec, rosmap_dec_adrelated, degree_limit = 1)
write_xlsx(ad_gene_network$edges, 'output/AD_related_network_edges.xlsx')

# 8. Protein Coding and lncRNA Analysis 
coding_gene_list <- read.delim("gene_with_protein_product.txt")
long_non_coding <- read.delim('RNA_long_non-coding.txt')
# Calculate intersections and perform Fisher's exact test here based on your original logic.
# Example:
# decreasing_table <- matrix(c(2972, 421, 12979, 9502), nrow = 2, byrow = TRUE)
# fisher_decreasing <- fisher.test(decreasing_table)
