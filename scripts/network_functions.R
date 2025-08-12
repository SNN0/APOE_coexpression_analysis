# File: network_functions.R
# This script contains a collection of functions for network analysis, including
# finding common gene pairs, exploring gene connections, and exporting results.

# 1. Load Libraries
library(data.table) # For efficient data manipulation
library(biomaRt)    # For gene ID conversion
library(writexl)    # For exporting data frames to Excel

# 2. Utility Functions

# Function to find common gene pairs between two datasets.
# It standardizes the pair order (e.g., alphabetically) before merging.
find_common_pairs <- function(data1, data2) {
  reorder_pairs <- function(data) {
    data[, .(gene1 = pmin(gene1, gene2), gene2 = pmax(gene1, gene2))]
  }
  
  ordered_data1 <- reorder_pairs(data1)
  ordered_data2 <- reorder_pairs(data2)
  
  common_pairs <- merge(ordered_data1, ordered_data2, by = c("gene1", "gene2"))
  
  return(list(common_pairs = common_pairs, count_common = nrow(common_pairs)))
}

# Function to convert gene symbols to Ensembl IDs.
convert_symbol_to_ensembl <- function(df, gene1_col, gene2_col) {
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene_symbols <- unique(c(df[[gene1_col]], df[[gene2_col]]))
  
  mapping <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id'),
                   filters = 'hgnc_symbol', 
                   values = gene_symbols, 
                   mart = ensembl)
  
  df <- merge(df, mapping, by.x = gene1_col, by.y = 'hgnc_symbol', all.x = TRUE)
  df <- merge(df, mapping, by.x = gene2_col, by.y = 'hgnc_symbol', all.x = TRUE, suffixes = c('_gene1', '_gene2'))
  
  return(df)
}

# Function to write a list of data frames to an Excel file.
write_all_dataframes_to_excel <- function(data_list, base_path, file_extension = ".xlsx") {
  if (!dir.exists(base_path)) {
    dir.create(base_path, recursive = TRUE)
    cat("Directory created:", base_path, "\n")
  }
  
  for (element_name in names(data_list)) {
    file_name <- paste0(element_name, file_extension)
    file_path <- file.path(base_path, file_name)
    writexl::write_xlsx(data_list[[element_name]], file_path)
    cat("Saved:", file_path, "\n")
  }
}

# 3. Network Exploration Functions

# Function to find all connections up to a specified degree for a single gene.
find_all_connections <- function(df, gene, degree_limit = Inf) {
  all_connections <- data.frame(gene1 = character(), gene2 = character(), stringsAsFactors = FALSE)
  current_genes <- c(gene)
  degree <- 1
  seen_genes <- c(gene)
  
  while (length(current_genes) > 0 & degree <= degree_limit) {
    new_connections <- df[df$gene1 %in% current_genes | df$gene2 %in% current_genes, ]
    all_connections <- unique(rbind(all_connections, new_connections))
    next_genes <- unique(c(new_connections$gene1, new_connections$gene2))
    current_genes <- setdiff(next_genes, seen_genes)
    seen_genes <- unique(c(seen_genes, current_genes))
    degree <- degree + 1
  }
  return(all_connections)
}

# Function to find all connections up to a specified degree for multiple genes.
# It returns a list of data frames, one for each seed gene.
find_all_connections_multiple_genes <- function(df, gene_list, degree_limit = Inf) {
  all_results <- list()
  
  for (gene in gene_list) {
    current_genes <- c(gene)
    seen_genes <- c(gene)
    degree <- 1
    
    while (length(current_genes) > 0 & degree <= degree_limit) {
      new_connections <- df[df$gene1 %in% current_genes | df$gene2 %in% current_genes, ]
      next_genes <- unique(c(new_connections$gene1, new_connections$gene2))
      current_genes <- setdiff(next_genes, seen_genes)
      seen_genes <- unique(c(seen_genes, current_genes))
      degree <- degree + 1
    }
    
    final_edges <- df[df$gene1 %in% seen_genes & df$gene2 %in% seen_genes, ]
    all_results[[gene]] <- final_edges
  }
  return(all_results)
}

# Function to find all connections for multiple genes and combine them into a single network.
find_all_connections_combined <- function(df, gene_list, degree_limit = Inf) {
  all_edges <- data.frame(gene1 = character(), gene2 = character(), stringsAsFactors = FALSE)
  all_nodes <- character()
  
  for (gene in gene_list) {
    current_genes <- c(gene)
    seen_genes <- c(gene)
    degree <- 1
    
    while (length(current_genes) > 0 && degree <= degree_limit) {
      new_connections <- df[df$gene1 %in% current_genes | df$gene2 %in% current_genes, ]
      all_edges <- unique(rbind(all_edges, new_connections))
      next_genes <- unique(c(new_connections$gene1, new_connections$gene2))
      current_genes <- setdiff(next_genes, seen_genes)
      seen_genes <- unique(c(seen_genes, current_genes))
      degree <- degree + 1
    }
    all_nodes <- union(all_nodes, seen_genes)
  }
  
  combined_edges <- df[df$gene1 %in% all_nodes & df$gene2 %in% all_nodes, ]
  all_nodes <- union(all_nodes, gene_list)
  gene_info <- data.frame(gene = all_nodes, stringsAsFactors = FALSE)
  gene_info$type <- ifelse(gene_info$gene %in% gene_list, "main", "side")
  
  return(list(edges = combined_edges, gene_info = gene_info))
}
