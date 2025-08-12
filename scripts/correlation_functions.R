# File: correlation_functions.R
# This script contains all the functions required for calculating correlation matrices
# and filtering co-expression patterns based on APOE genotypes.
# The functions are structured to precisely match the user's original logic.

# 1. Load Libraries
library(Hmisc)
library(data.table)
library(foreach)
library(doParallel)
library(dplyr)

# 2. Helper Functions

# Helper function to perform BH-adjusted p-value calculation.
adjust_p_values_bh <- function(p_matrix) {
  p_values_vector <- as.vector(p_matrix)
  p_values_vector <- p_values_vector[!is.na(p_values_vector)]
  adjusted_p_values <- p.adjust(p_values_vector, method = 'BH')
  
  adjusted_p_matrix <- p_matrix
  adjusted_p_matrix[!is.na(adjusted_p_matrix)] <- adjusted_p_values
  return(adjusted_p_matrix)
}

# 3. Main Correlation Functions

# Function 1: Calculate correlation matrices for each APOE genotype using parallel processing.
# It uses the VST-transformed expression data and filters by genotype from the metadata.
calculate_correlations <- function(expression_data, meta_data, genotypes) {
  correlations <- list()
  
  # Setup parallel processing
  n_cores <- 4  # Adjust as needed, but do not exceed available cores
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Use foreach to run correlation calculation for each genotype in parallel
  results <- foreach(genotype = genotypes, .packages = c('Hmisc', 'data.table', 'dplyr')) %dopar% {
    subset_data <- expression_data[meta_data$apoe_genotype == genotype, ]
    corr_result <- Hmisc::rcorr(as.matrix(subset_data), type = 'pearson')
    adj_p_matrix <- adjust_p_values_bh(corr_result$P)
    list(correlations = corr_result$r, p_value = adj_p_matrix)
  }
  
  stopCluster(cl)
  names(results) <- genotypes
  return(results)
}

# Function 2: Filter significant correlation pairs based on correlation and p-value thresholds.
# This function directly implements the user's `filter_signi_corr_optimized` logic.
filter_significant_pairs <- function(correlations_results, genotype, threshold_corr = 0.95, threshold_p = 0.01) {
  ind <- which(upper.tri(correlations_results[[genotype]]$correlations, diag = FALSE), arr.ind = TRUE)
  
  signi_pairs_dt <- data.table(
    gene1 = dimnames(correlations_results[[genotype]]$correlations)[[2]][ind[, 2]],
    gene2 = dimnames(correlations_results[[genotype]]$correlations)[[1]][ind[, 1]],
    corr_val = correlations_results[[genotype]]$correlations[ind],
    p_val = correlations_results[[genotype]]$p_value[ind]
  )
  
  filtered_pairs <- signi_pairs_dt[abs(corr_val) > threshold_corr & p_val < threshold_p]
  return(setNames(list(filtered_pairs), genotype))
}

# Function 3: Identify decreasing and increasing patterns, referencing E2E3 as the base.
# This function directly implements the user's `filter_correlation_pattern_all_optimized` logic.
filter_correlation_patterns_E2E3_ref <- function(significant_pairs, correlations, genotypes, threshold_diff = 0.15) {
  ref_genotype <- genotypes[1]  # E2E3 as reference
  mid_genotype <- genotypes[2]  # E3E3
  target_genotype <- genotypes[3]  # E3E4
  
  sig_pairs_dt <- data.table(
    gene1 = significant_pairs[[ref_genotype]]$gene1,
    gene2 = significant_pairs[[ref_genotype]]$gene2
  )
  
  sig_pairs_dt[, ref_corr := correlations[[ref_genotype]]$correlations[cbind(gene1, gene2)]]
  sig_pairs_dt[, mid_corr := correlations[[mid_genotype]]$correlations[cbind(gene1, gene2)]]
  sig_pairs_dt[, target_corr := correlations[[target_genotype]]$correlations[cbind(gene1, gene2)]]
  
  sig_pairs_dt[, diff_corr_1 := ref_corr - mid_corr]
  sig_pairs_dt[, diff_corr_2 := mid_corr - target_corr]
  
  # Decreasing pattern logic
  sig_pairs_dt[, pass_decreasing := (ref_corr > 0 & diff_corr_1 >= threshold_diff & diff_corr_2 >= threshold_diff & mid_corr > 0 & target_corr > 0) |
                 (ref_corr < 0 & diff_corr_1 <= -threshold_diff & diff_corr_2 <= -threshold_diff & mid_corr < 0 & target_corr < 0)]
  
  # Increasing pattern logic
  sig_pairs_dt[, pass_increasing := (ref_corr > 0 & diff_corr_1 <= -threshold_diff & diff_corr_2 <= -threshold_diff & mid_corr > 0 & target_corr > 0) |
                 (ref_corr < 0 & diff_corr_1 >= threshold_diff & diff_corr_2 >= threshold_diff & mid_corr < 0 & target_corr < 0)]
  
  positive_cor_decreasing <- sig_pairs_dt[pass_decreasing == TRUE & ref_corr > 0, .(gene1, gene2, ref_corr)]
  negative_cor_decreasing <- sig_pairs_dt[pass_decreasing == TRUE & ref_corr < 0, .(gene1, gene2, ref_corr)]
  positive_cor_increasing <- sig_pairs_dt[pass_increasing == TRUE & ref_corr > 0, .(gene1, gene2, ref_corr)]
  negative_cor_increasing <- sig_pairs_dt[pass_increasing == TRUE & ref_corr < 0, .(gene1, gene2, ref_corr)]
  
  return(list(
    positive_cor_decreasing = positive_cor_decreasing,
    negative_cor_decreasing = negative_cor_decreasing,
    positive_cor_increasing = positive_cor_increasing,
    negative_cor_increasing = negative_cor_increasing
  ))
}

# Function 4: Identify increasing patterns, referencing E3E4 as the base.
# This function directly implements the user's `filter_correlation_pattern_increasing_E3E4` logic.
filter_correlation_patterns_E3E4_ref <- function(significant_pairs, correlations, genotypes, threshold_diff = 0.15) {
  ref_genotype <- genotypes[3]  # E3E4 as reference
  mid_genotype <- genotypes[2]  # E3E3
  target_genotype <- genotypes[1]  # E2E3
  
  sig_pairs_dt <- data.table(
    gene1 = significant_pairs[[ref_genotype]]$gene1,
    gene2 = significant_pairs[[ref_genotype]]$gene2
  )
  
  sig_pairs_dt[, ref_corr := correlations[[ref_genotype]]$correlations[cbind(gene1, gene2)]]
  sig_pairs_dt[, mid_corr := correlations[[mid_genotype]]$correlations[cbind(gene1, gene2)]]
  sig_pairs_dt[, target_corr := correlations[[target_genotype]]$correlations[cbind(gene1, gene2)]]
  
  sig_pairs_dt[, diff_corr_1 := ref_corr - mid_corr]
  sig_pairs_dt[, diff_corr_2 := mid_corr - target_corr]
  
  # Increasing pattern logic (without crossing zero)
  sig_pairs_dt[, pass_increasing := (ref_corr > 0 & diff_corr_1 >= threshold_diff & diff_corr_2 >= threshold_diff & mid_corr > 0 & target_corr > 0) |
                 (ref_corr < 0 & diff_corr_1 <= -threshold_diff & diff_corr_2 <= -threshold_diff & mid_corr < 0 & target_corr < 0)]
  
  positive_cor_increasing <- sig_pairs_dt[pass_increasing == TRUE & ref_corr > 0, .(gene1, gene2, ref_corr)]
  negative_cor_increasing <- sig_pairs_dt[pass_increasing == TRUE & ref_corr < 0, .(gene1, gene2, ref_corr)]
  
  return(list(
    positive_cor_increasing = positive_cor_increasing,
    negative_cor_increasing = negative_cor_increasing
  ))
}

# 4. Visualization Functions
# 
gene_mapping <- c(
  "ENSG00000130203" = "APOE",
  "ENSG00000176340" = "COX8A",
  "ENSG00000087258" = "GNAO1",
  "ENSG00000111237" = "VPS29",
  "ENSG00000166471" = "TMEM41B",
  "ENSG00000139597" = "N4BP2L1"
)

# Fonksiyon tan??mlama
plot_gene_correlation <- function(meta_data, expression_data, gene1, gene2) {
  
  # Meta veriyi specimenID'ye g??re e??le??tir
  meta_data$specimenID <- as.character(meta_data$specimenID)
  
  # Se??ilen iki genin ekspresyon verisini al
  expression_values <- expression_data %>%
    rownames_to_column(var = "specimenID") %>%
    select(specimenID, all_of(c(gene1, gene2)))
  
  # Meta veriye ekspresyon verisini ekle
  merged_data <- inner_join(meta_data, expression_values, by = "specimenID")
  
  # E??imi hesaplamak i??in lineer modelleri olu??tur
  slope_data <- merged_data %>%
    group_by(apoe_genotype) %>%
    summarize(slope = coef(lm(.data[[gene2]] ~ .data[[gene1]]))[2])  # Regresyon katsay??s??n?? al
  
  # **Gen isimlerini al 
  gene1_name <- ifelse(gene1 %in% names(gene_mapping), gene_mapping[[gene1]], gene1)
  gene2_name <- ifelse(gene2 %in% names(gene_mapping), gene_mapping[[gene2]], gene2)
  
  # Facet grid scatter plot ??izimi
  ggplot(merged_data, aes(x = .data[[gene1]], y = .data[[gene2]], color = apoe_genotype)) +
    geom_point(aes(shape = apoe_genotype), size = 4, alpha = 0.8, stroke = 1) +  # B??y??k ve belirgin noktalar
    geom_smooth(method = "lm", se = TRUE, linetype = "solid", linewidth = 1, alpha = 0.6) +  # Regresyon ??izgisi
    facet_wrap(~apoe_genotype) +  # Genotipe g??re ay??r
    scale_color_manual(values = c("E2E3" = "#1f77b4", "E3E3" = "#ff7f0e", "E3E4" = "#2ca02c")) +  # ??zel renk paleti
    scale_shape_manual(values = c("E2E3" = 16, "E3E3" = 17, "E3E4" = 15)) +  # Farkl?? ??ekiller
    theme_minimal(base_size = 14) +  # Modern ve temiz bir tema
    theme(
      panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 14, face = "bold"),
      legend.position = "top"
    ) +
    labs(
      title = paste(gene1_name, "vs", gene2_name, "Expression by APOE Genotype"),
      subtitle = "Slope (regression coefficient) indicates correlation strength",
      x = paste(gene1_name, "Expression"),
      y = paste(gene2_name, "Expression"),
      color = "APOE Genotype",
      shape = "APOE Genotype"
    ) +
    geom_text(data = slope_data, aes(x = Inf, y = -Inf, label = paste0("Slope: ", round(slope, 3))),
              hjust = 1.1, vjust = -1.1, size = 5, color = "black", inherit.aes = FALSE)
}

plot_gene_correlation(meta.data, as.data.frame(norm.count.vst), "ENSG00000130203", "ENSG00000087258")
