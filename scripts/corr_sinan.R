
data=readRDS('main_data.rds')

count.data = data$count_data
meta.data = data$meta_data

library(DESeq2)
library(tidyverse)

dds = DESeqDataSetFromMatrix(countData = count.data,
                             colData = meta.data,
                             design = ~apoe_genotype + age_death +msex +pmi)

dds = DESeq(dds)

resultsNames(dds)

norm.count = counts(dds, normalized = T)


dds.vst = vst(dds) #log2 
norm.count.vst = assay(dds.vst) %>% t()

#design matrix
design.matrix = model.matrix(~+msex+pmi+age_death, meta.data)

long_non_coding = read.delim('RNA_long_non-coding.txt', header = T, sep = '\t',stringsAsFactors = F)
long_non_coding = long_non_coding$ensembl_gene_id

sum(rownames(count.data) %in% long_non_coding)

# norm.count.vst'yi data frame'e ??evir ve sat??r adlar??n?? sample ID'leri olarak ayarla
norm_data11<- as.data.frame(norm.count.vst)
rownames(norm_data11) <- colnames(norm.count.vst)


# Verileri birle??tir 
combined_data <- merge(norm_data11, meta.data, by.x = "row.names", by.y = "specimenID")
rownames(combined_data) = combined_data$Row.names

gen_expression= combined_data[, c("ENSG00000000003","ENSG00000000419", "apoe_genotype")]


rep(gen_expression$apoe_genotype, 2)

meta.data$specimenID <- as.character(meta.data$specimenID)

# Gerekli k??t??phaneleri y??kle
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggthemes)

# Ensembl ID - Gen ??smi E??leme
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



# Fonksiyon tan??mlama
plot_two_genes_boxplot <- function(meta_data, expression_data, gene1, gene2) {
  
  # Meta veriyi specimenID'ye g??re e??le??tir
  meta_data$specimenID <- as.character(meta_data$specimenID)
  
  # Se??ilen iki genin ekspresyon verisini al
  expression_values <- expression_data %>%
    rownames_to_column(var = "specimenID") %>%
    select(specimenID, all_of(c(gene1, gene2)))
  
  # Meta veriye ekspresyon verisini ekle
  merged_data <- inner_join(meta_data, expression_values, by = "specimenID")
  
  # Uzun format (long format) d??n??????m?? (ggplot i??in gerekli)
  long_data <- merged_data %>%
    pivot_longer(cols = c(gene1, gene2), names_to = "Gene", values_to = "Expression")
  
  # Boxplot ??izimi
  ggplot(long_data, aes(x = apoe_genotype, y = Expression, fill = Gene)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Boxplot (Medyan ve IQR)
    geom_jitter(aes(color = Gene), shape = 21, size = 2, width = 0.2, alpha = 0.7) +  # Hafif kayd??r??lm???? noktalar
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, aes(fill = Gene), color = "black") +  # Ortalama noktas??
    theme_minimal() +
    labs(title = paste(gene1, "vs", gene2, "Expression by APOE Genotype"),
         x = "APOE Genotype",
         y = "Expression Level") +
    scale_fill_manual(values = c(gene1 = "blue", gene2 = "red")) +
    scale_color_manual(values = c(gene1 = "blue", gene2 = "red")) +
    theme(legend.position = "top")
}

# ??rnek kullan??m
plot_two_genes_boxplot(meta.data, as.data.frame(norm.count.vst), "ENSG00000130203", "ENSG00000139597")


#---------Optimize---------------


library(Hmisc)


cor_with_pvalues_yeni = function(data, method) {
  
  # rcorr fonksiyonunu kullanarak korelasyonlar?? ve p-de??erlerini hesaplama
  result = Hmisc::rcorr(as.matrix(data), type = method)
  
  # Korelasyon ve p-de??erleri
  corr_matrix = result$r
  p_matrix = result$P
  
  # BH (Benjamini-Hochberg) ile d??zeltilmi?? p-de??erlerini hesaplama
  # rcorr ????kt??s??ndaki p-de??erlerini vekt??re ??eviriyoruz ve NA de??erleri ????kar??yoruz
  p_values_vector = as.vector(p_matrix)
  p_values_vector = p_values_vector[!is.na(p_values_vector)]
  
  # BH y??ntemiyle p-de??erlerini d??zeltiyoruz
  adjusted_p_values = p.adjust(p_values_vector, method = "BH")
  
  # D??zeltilmi?? p-de??erlerini yeni bir p_matrix olu??turuyoruz
  adjusted_p_matrix = p_matrix  # Ayn?? boyutta bir matris olu??turuyoruz
  adjusted_p_matrix[!is.na(adjusted_p_matrix)] = adjusted_p_values
  
  # Sonu??lar?? liste halinde d??nd??rmek
  list(correlations = corr_matrix, p_value = adjusted_p_matrix)
  
}


#----------------------------------

cor_with_pvalues = function(data) {
  n = ncol(data)
  corr_matrix = cor(data, method = "pearson")
  p_matrix = matrix(NA, n, n)
  diag(p_matrix) = 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      test_result = cor.test(data[, i], data[, j], method = "pearson")
      p_matrix[i, j] = p_matrix[j, i] <- test_result$p.value
    }
  }
  list(correlations = corr_matrix, p_values = p_matrix)
}

library(foreach)
library(doParallel)

cal_corr = function(data, genotype) {
  correlations = list()
  
  results = foreach(genotype = genotypes, .packages = c('foreach', 'doParallel'), .export = c('meta.data', 'cor_with_pvalues_yeni')) %dopar%  {
    subset_data = data[meta.data$apoe_genotype == genotype, ]
    corr_matrix = cor_with_pvalues_yeni(subset_data, 'pearson')
    correlations[[genotype]] = corr_matrix
  }
  
  names(results) = genotypes
  return(results)
  
}

#n_cores = 2 #max 32
#cl = makeCluster(n_cores)
#registerDoParallel(cl)

genotypes = c('E2E3', 'E3E3', 'E3E4')

#correlations = cal_corr(norm.count.vst, genotypes)

#stopCluster(cl)

filter_signi_corr = function(correlations_results, threshold_corr = 0.95, threshold_p = 0.01) {
  
  signi_pairs = list()
  
  for(genotype in names(correlations_results)){
    
    ind = which(upper.tri(correlations_results[[genotype]]$correlations, diag=F), arr.ind = TRUE)
    bc=data.frame( col = dimnames(correlations_results[[genotype]]$correlations)[[2]][ind[,2]] , # ind 2 col numaralar??
                   row = dimnames(correlations_results[[genotype]]$correlations)[[1]][ind[,1]] , # ind 1 row
                   corr_val = correlations_results[[genotype]]$correlations[ind],
                   p_val = correlations_results[[genotype]]$p_values[ind])
    
    bc = bc %>% filter(., abs(corr_val) > threshold_corr)
    bc = bc %>% filter(., p_val < threshold_p)
    
    signi_pairs[[genotype]] = bc
    
  }
  
  return(signi_pairs)
}

#filter_signi_cor E2E3

#test edildi.
filter_signi_corr_E2E3 = function(correlations_results, threshold_corr = 0.95, threshold_p = 0.01) {
  
  signi_pairs = list()
  genotype = 'E2E3'
  
  
    
  ind = which(upper.tri(correlations_results[[genotype]]$correlations, diag=F), arr.ind = TRUE)
 
  bc=data.frame( col = dimnames(correlations_results[[genotype]]$correlations)[[2]][ind[,2]] , # ind 2 col numaralar??
                   row = dimnames(correlations_results[[genotype]]$correlations)[[1]][ind[,1]] , # ind 1 row
                   corr_val = correlations_results[[genotype]]$correlations[ind],
                   p_val = correlations_results[[genotype]]$adjusted_p_values[ind])
    
  bc = bc %>% filter(., abs(corr_val) > threshold_corr)
  bc = bc %>% filter(., p_val < threshold_p)
    
  signi_pairs[[genotype]] = bc
    
  
  
  return(signi_pairs)
}



#------------------------------------------------------------------------------------------------------------------

library(data.table)
#optimize edilmi??
filter_signi_corr_optimized = function(correlations_results, genotype, threshold_corr = 0.95, threshold_p = 0.01) {
  
  # Extract significant correlation pairs using the upper triangle
  ind = which(upper.tri(correlations_results[[genotype]]$correlations, diag=F), arr.ind = TRUE)
  
  # Convert the result into a data.table for faster processing
  bc = data.table(
    gene1 = dimnames(correlations_results[[genotype]]$correlations)[[2]][ind[,2]],
    gene2 = dimnames(correlations_results[[genotype]]$correlations)[[1]][ind[,1]],
    corr_val = correlations_results[[genotype]]$correlations[ind],
    p_val = correlations_results[[genotype]]$p_value[ind]
  )
  
  # Filter based on the correlation and p-value thresholds
  bc = bc[abs(corr_val) > threshold_corr & p_val < threshold_p]
  
  return(setNames(list(bc), genotype))
}


#filter pattern optimize ??al????t??.
library(data.table)

filter_correlation_pattern_all_optimized = function(significant_pairs, correlations, threshold_diff = 0.15) {
  
  ref_genotype = genotypes[1]  # E2E3 as reference
  mid_genotype = genotypes[2]  # E3E3 as middle group
  target_genotype = genotypes[3]  # E3E4 as the last group
  
  # Extract significant pairs into a data.table
  sig_pairs_dt = data.table(
    gene1 = significant_pairs[[ref_genotype]]$gene1,
    gene2 = significant_pairs[[ref_genotype]]$gene2
  )
  
  # Retrieve correlation values for all pairs at once using vectorized operations
  sig_pairs_dt[, ref_corr := correlations[[ref_genotype]]$correlations[cbind(gene1, gene2)]]
  sig_pairs_dt[, mid_corr := correlations[[mid_genotype]]$correlations[cbind(gene1, gene2)]]
  sig_pairs_dt[, target_corr := correlations[[target_genotype]]$correlations[cbind(gene1, gene2)]]
  
  # Calculate the differences between correlations
  sig_pairs_dt[, diff_corr_1 := ref_corr - mid_corr]
  sig_pairs_dt[, diff_corr_2 := mid_corr - target_corr]
  
  # Identify decreasing patterns
  sig_pairs_dt[, pass_decreasing := (ref_corr > 0 & diff_corr_1 >= threshold_diff & diff_corr_2 >= threshold_diff & mid_corr >0 & target_corr > 0 ) |
                 (ref_corr < 0 & diff_corr_1 <= -threshold_diff & diff_corr_2 <= -threshold_diff & mid_corr <0 & target_corr < 0)]
  
  # Identify increasing patterns
  sig_pairs_dt[, pass_increasing := (ref_corr > 0 & diff_corr_1 <= -threshold_diff & diff_corr_2 <= -threshold_diff & mid_corr >0 & target_corr > 0) |
                 (ref_corr < 0 & diff_corr_1 >= threshold_diff & diff_corr_2 >= threshold_diff & mid_corr <0 & target_corr < 0)]
  
  # Split the results into separate lists
  positive_cor_decreasing = sig_pairs_dt[ref_corr > 0 & pass_decreasing == TRUE, .(gene1, gene2, ref_corr)]
  positive_cor_increasing = sig_pairs_dt[ref_corr > 0 & pass_increasing == TRUE, .(gene1, gene2, ref_corr)]
  negative_cor_decreasing = sig_pairs_dt[ref_corr < 0 & pass_decreasing == TRUE, .(gene1, gene2, ref_corr)]
  negative_cor_increasing = sig_pairs_dt[ref_corr < 0 & pass_increasing == TRUE, .(gene1, gene2, ref_corr)]
  
  return(list(
    positive_cor_decreasing = positive_cor_decreasing,
    positive_cor_increasing = positive_cor_increasing,
    negative_cor_decreasing = negative_cor_decreasing,
    negative_cor_increasing = negative_cor_increasing
  ))
}

filter_correlation_pattern_increasing_E3E4 = function(significant_pairs, correlations, threshold_diff = 0.15) {
  
  ref_genotype = genotypes[3]  # E3E4 as reference
  mid_genotype = genotypes[2]  # E3E3 as middle group
  target_genotype = genotypes[1]  # E2E3 as the first group
  
  # Extract significant pairs into a data.table
  sig_pairs_dt = data.table(
    gene1 = significant_pairs[[ref_genotype]]$gene1,
    gene2 = significant_pairs[[ref_genotype]]$gene2
  )
  
  # Retrieve correlation values for all pairs at once using vectorized operations
  sig_pairs_dt[, ref_corr := correlations[[ref_genotype]]$correlations[cbind(gene1, gene2)]]
  sig_pairs_dt[, mid_corr := correlations[[mid_genotype]]$correlations[cbind(gene1, gene2)]]
  sig_pairs_dt[, target_corr := correlations[[target_genotype]]$correlations[cbind(gene1, gene2)]]
  
  # Calculate the differences between correlations (decreasing pattern from E3E4 to E2E3)
  sig_pairs_dt[, diff_corr_1 := ref_corr - mid_corr]
  sig_pairs_dt[, diff_corr_2 := mid_corr - target_corr]
  
  # Identify increasing patterns (without crossing zero)
  sig_pairs_dt[, pass_increasing := (ref_corr > 0 & diff_corr_1 >= threshold_diff & diff_corr_2 >= threshold_diff & mid_corr >0 & target_corr > 0 ) |
                 (ref_corr < 0 & diff_corr_1 <= -threshold_diff & diff_corr_2 <= -threshold_diff & mid_corr <0 & target_corr < 0)]
  
  # Split the results into increasing patterns
  positive_cor_increasing = sig_pairs_dt[ref_corr > 0 & pass_increasing == TRUE, .(gene1, gene2, ref_corr)]
  negative_cor_increasing = sig_pairs_dt[ref_corr < 0 & pass_increasing == TRUE, .(gene1, gene2, ref_corr)]
  
  return(list(
    positive_cor_increasing = positive_cor_increasing,
    negative_cor_increasing = negative_cor_increasing
  ))
}


#pval 0.01 ve 0.05 te  0.25 arat??kl?? e2e3 ten artan yok adj pval=0.05 te en kucuk |cor| = 0.818


signi_pairs_rosmap_e2e3 = filter_signi_corr_optimized(correlations, 'E2E3', threshold_corr = 0.95, threshold_p = 0.05)
pattern_rosmap_0.05_dec = filter_correlation_pattern_all_optimized(signi_pairs_rosmap_e2e3, correlations, threshold_diff = 0.25)




signi_pairs_rosmap_e3e4 = filter_signi_corr_optimized(correlation, 'E3E4', threshold_corr = 0.80, threshold_p = 0.05)
pattern_rosmap_0.05_inc = filter_correlation_pattern_increasing_E3E4(signi_pairs_rosmap_e3e4, correlation, threshold_diff = 0.25)



saveRDS(pattern_rosmap_0.05_dec, 'pattern/Rosmap_0.90cor_decreasing.rds')
saveRDS(pattern_rosmap_0.05_dec, 'pattern/Rosmap_0.95cor_decreasing.rds')

saveRDS(pattern_rosmap_0.05_inc, 'pattern/Rosmap_0.80cor_increasing.rds')




#inc - 0.9 cor,  0.25 dif,  => 0pos , 0neg padj 0.05

#inc - 0.9 cor,  0.2 dif,  => 0pos , 0neg padj 0.05

#inc - 0.85 cor,  0.25 dif,  => 28pos , 18neg padj 0.05

#inc - 0.85 cor,  0.2 dif,  => 94pos , 46neg padj 0.05

#inc - 0.8 cor,  0.25 dif,  => 341pos , 215neg padj 0.05

#inc - 0.8 cor,  0.2 dif,  => 1077 pos , 705 neg  padj 0.05

#inc - 0.75 cor,  0.25 dif,  => 1556pos , 1371neg padj 0.05

#inc - 0.75 cor,  0.2 dif,  => 4905 pos , 4120 neg  padj 0.05


#signi E2E3 0.05

#signi olanlar?? filtrelemek i??in
signi_pairs = filter_signi_corr(correlations, threshold_p = 0.01)

#artan i??in refi 0.05 cor olarak ayarlard??m ... 
signi_pairs_0.05 = filter_signi_corr(correlations,threshold_corr = 0.05,threshold_p = 0.01)

# 0.04 te 10 tane ????kmas?? laz??m kontrol et. pos increasing de 
filter_correlation_pattern_all = filter_correlation_pattern_all(signi_pairs_0.05, correlations, threshold_diff = 0.04)




signi_pairs_normal = filter_signi_corr(correlations, threshold_corr=0.95, threshold_p = 0.01)
signi_pairs_increasing = filter_signi_corr(correlations, threshold_corr = 0, threshold_p= 0.01)

signi_pairs_pval0.05AGrubu = filter_signi_corr(correlations, threshold_corr = 0, threshold_p = 0.05)




#deneme en yenisi :
filter_pattern_deneme = filter_correlation_pattern_all_yeni(signi_pairs_normal,signi_pairs_increasing, correlations, threshold_diff =0.25)


#A grubundan artan var m??? cor 0
filter_pattern_Agrup_increasing = filter_correlation_pattern_all(signi_pairs_pval0.05AGrubu, correlations, threshold_diff = 0.25)





#en yeni yaz??lan fonksiyon
filter_correlation_pattern_all_yeni = function(significant_pairs,signi_increasing, correlations, threshold_diff = 0.25) {
  positive_cor_decreasing = list()
  positive_cor_increasing = list()
  negative_cor_decreasing = list()
  negative_cor_increasing = list()
  
  ref_genotype = genotypes[1]  # E2E3 as reference for decreasing pattern
  target_genotype = genotypes[3]  # E3E4 as reference for increasing pattern
  
  # Decreasing pattern check from E2E3 -> E3E3 -> E3E4
  for (pair in 1:nrow(significant_pairs[[ref_genotype]])) {
    gene1 = significant_pairs[[ref_genotype]][pair,2]
    gene2 = significant_pairs[[ref_genotype]][pair,1]
    
    ref_corr = correlations[[ref_genotype]]$correlations[gene1, gene2]
    mid_corr = correlations[[genotypes[2]]]$correlations[gene1, gene2]
    target_corr = correlations[[target_genotype]]$correlations[gene1, gene2]
    
    diff_corr_1 = ref_corr - mid_corr
    diff_corr_2 = mid_corr - target_corr
    
    pass_pattern_decreasing = TRUE
    
    if (ref_corr > 0) {
      if (diff_corr_1 < threshold_diff || diff_corr_2 < threshold_diff) {
        pass_pattern_decreasing = FALSE
      }
    } else {
      if (diff_corr_1 > -threshold_diff || diff_corr_2 > -threshold_diff) {
        pass_pattern_decreasing = FALSE
      }
    }
    
    if (ref_corr > 0 && pass_pattern_decreasing) {
      positive_cor_decreasing[[paste(gene1, gene2, sep = "_")]] = ref_corr
    } else if (ref_corr < 0 && pass_pattern_decreasing) {
      negative_cor_decreasing[[paste(gene1, gene2, sep = "_")]] = ref_corr
    }
  }
  
  # Increasing pattern check from E3E4 -> E3E3 -> E2E3
  for (pair in 1:nrow(signi_increasing[[target_genotype]])) {
    gene1 = signi_increasing[[target_genotype]][pair,2]
    gene2 = signi_increasing[[target_genotype]][pair,1]
    
    ref_corr = correlations[[ref_genotype]]$correlations[gene1, gene2]
    mid_corr = correlations[[genotypes[2]]]$correlations[gene1, gene2]
    target_corr = correlations[[target_genotype]]$correlations[gene1, gene2]
    
    diff_corr_1 = target_corr - mid_corr
    diff_corr_2 = mid_corr - ref_corr
    
    pass_pattern_increasing = TRUE
    
    if (target_corr > 0) {
      if (diff_corr_1 < threshold_diff || diff_corr_2 < threshold_diff || ref_corr <= 0) {
        pass_pattern_increasing = FALSE
      }
    } else {
      if (diff_corr_1 > -threshold_diff || diff_corr_2 > -threshold_diff || ref_corr >= 0) {
        pass_pattern_increasing = FALSE
      }
    }
    
    if (target_corr > 0 && pass_pattern_increasing) {
      positive_cor_increasing[[paste(gene1, gene2, sep = "_")]] = target_corr
    } else if (target_corr < 0 && pass_pattern_increasing) {
      negative_cor_increasing[[paste(gene1, gene2, sep = "_")]] = target_corr
    }
  }
  
  return(list(
    positive_cor_decreasing = positive_cor_decreasing,
    positive_cor_increasing = positive_cor_increasing,
    negative_cor_decreasing = negative_cor_decreasing,
    negative_cor_increasing = negative_cor_increasing
  ))
}




## DENEME **
#daha yeni pattern: kontrol edildi.Sonu??lar??n kulln??ld??g?? fonksiyon

filter_correlation_pattern_all = function(significant_pairs, correlations, threshold_diff = 0.15) {
  positive_cor_decreasing = list()
  positive_cor_increasing = list()
  negative_cor_decreasing = list()
  negative_cor_increasing = list()
  
  ref_genotype = genotypes[1]  # E2E3 as reference
  mid_genotype = genotypes[2]  # E3E3 as middle group
  target_genotype = genotypes[3]  # E3E4 as the last group
  
  for (pair in 1:nrow(significant_pairs[[ref_genotype]])) {
    gene1 = significant_pairs[[ref_genotype]][pair,2]
    gene2 = significant_pairs[[ref_genotype]][pair,1]
    
    ref_corr = correlations[[ref_genotype]]$correlations[gene1, gene2]
    mid_corr = correlations[[mid_genotype]]$correlations[gene1, gene2]
    target_corr = correlations[[target_genotype]]$correlations[gene1, gene2]
    
    diff_corr_1 = ref_corr - mid_corr
    diff_corr_2 = mid_corr - target_corr
    
    pass_pattern_decreasing = TRUE
    pass_pattern_increasing = TRUE
    
    if (ref_corr > 0) {
      # Check for decreasing pattern
      if (diff_corr_1 < threshold_diff || diff_corr_2 < threshold_diff) {
        pass_pattern_decreasing = FALSE
      }
      # Check for increasing pattern
      if (diff_corr_1 > -threshold_diff  || diff_corr_2 > -threshold_diff  ) { # 
        pass_pattern_increasing = FALSE
      }
    } else {
      # Check for decreasing pattern
      if (diff_corr_1 > -threshold_diff || diff_corr_2 > -threshold_diff) {
        pass_pattern_decreasing = FALSE
      }
      # Check for increasing pattern
      if (diff_corr_1 < threshold_diff || diff_corr_2 < threshold_diff) { # 
        pass_pattern_increasing = FALSE
      }
    }
    
    if (ref_corr > 0) {
      if (pass_pattern_decreasing) {
        positive_cor_decreasing[[paste(gene1, gene2, sep = "_")]] = ref_corr
      }
      if (pass_pattern_increasing) {
        positive_cor_increasing[[paste(gene1, gene2, sep = "_")]] = ref_corr
      }
    } else {
      if (pass_pattern_decreasing) {
        negative_cor_decreasing[[paste(gene1, gene2, sep = "_")]] = ref_corr
      }
      if (pass_pattern_increasing) {
        negative_cor_increasing[[paste(gene1, gene2, sep = "_")]] = ref_corr
      }
    }
  }
  
  return(list(
    positive_cor_decreasing = positive_cor_decreasing,
    positive_cor_increasing = positive_cor_increasing,
    negative_cor_decreasing = negative_cor_decreasing,
    negative_cor_increasing = negative_cor_increasing
  ))
}












filtered_patterns = filter_correlation_pattern101(signi_pairs,correlations,0.25) # 3926 pos_dec,  3360 neg_dec
filtered_patterns0.35 = filter_correlation_pattern101(signi_pairs,correlations,0.35) #pos_dec 1303, 1067 neg_dec


# 
# p-val matrixe row ve column isimleri yok... 
for(genotype in genotypes) {
  p_matrix = correlations[[genotype]]$p_values
  colnames(p_matrix) = colnames(correlations$E2E3$correlations)
  rownames(p_matrix) = colnames(correlations$E2E3$correlations)
  correlations[[genotype]]$p_values = p_matrix
}

#adj- pval - bu i??lemi yapt??n, correlation.rds  p.adjli
adjust_p_values = function(p_matrix){
  p_vector = as.vector(p_matrix)
  adjusted_p_vector = p.adjust(p_vector, method = 'BH')
  matrix(adjusted_p_vector, nrow = nrow(p_matrix), ncol = ncol(p_matrix))
}

correlations = lapply(correlations, function(result){
  result$p_values = adjust_p_values(result$p_values)
  result
})


#### -------------------------------

select_top_pairs = function(correlations_pattern, correlations_results, top_n = 5){
  
  top_pairs = list()
  
  for(pattern_name in names(correlations_pattern)){
    pattern= correlations_pattern[[pattern_name]]
    
    if(length(pattern) > 0){ #pattern bo?? mu ? de??ilse devam
      pattern_df = data.frame(
        gene_pair = names(pattern),
        correlation = unlist(pattern)
      )
      
      #p-val ekle
      pattern_df$p_value = sapply(pattern_df$gene_pair, function(pair){
        
        genes = unlist(strsplit(pair, '_')) #gen ??ift isimlerini ay??r
        correlations_results[['E2E3']]$p_values[genes[1],genes[2]] # pval ??ek
        
      })
      # en y??ksek cor ve en d??????k p-val
      pattern_df = pattern_df[order(-pattern_df$correlation, pattern_df$p_value),]
      
      top_pairs[[pattern_name]] = head(pattern_df, top_n)
    }
  }
  
  return(top_pairs)
}

top_pairs=select_top_pairs(filtered_patterns, correlations,5)

#belirli bir gen ismim ile ??ift bulma

find_gene_pairs = function(gene_name, filtered_patterns){
  gene_pairs = list()
  
  for(pattern_name in names(filtered_patterns)){
    pattern = filtered_patterns[[pattern_name]]
    
    #belirli geni i??eren geni i??eren ??ift
    matching_pairs = pattern[grep(gene_name, names(pattern))]
    gene_pairs[[pattern_name]] = matching_pairs
  }
  
  return(gene_pairs)
  
  
}

find_gene_pairs('ENSG00000130203',filteredPattern)
# APOE nin oldu??u ??iftler:
# pos_cor_dec :  ENSG00000087258, ENSG00000176340
# neg cor dec : ENSG00000111237, ENSG00000139597, ENSG00000145388, ENSG00000166471

find_gene_pairs('ENSG00000002330',filteredPattern) #BAD gene
# pos BAD: ENSG00000061337, ENSG00000161016,  ENSG00000203950, ENSG00000243279
#neg BAD:  ENSG00000133114

find_gene_pairs('ENSG00000176340', filteredPattern) #COX8A



COX8A_pos_pair = names(find_gene_pairs('ENSG00000176340', filteredPattern)$positive_cor_decreasing)
COX8A_neg_pair = names(find_gene_pairs('ENSG00000176340', filteredPattern)$negative_cor_decreasing)


APOE_pos_pair = c('ENSG00000087258_ENSG00000130203', 'ENSG00000130203_ENSG00000176340')
APOE_neg_pair = c('ENSG00000111237_ENSG00000130203', 'ENSG00000130203_ENSG00000139597', 'ENSG00000130203_ENSG00000166471')


plot_selected_pair = function(selected_pairs, correlations){
  #isted??im gen ??ifti grafi??ini ??izme liste olarak vericem
  
  data_list = lapply(selected_pairs, function(gene_pair){
    
    genes = unlist(strsplit(gene_pair, '_'))
    gene1 = genes[1]
    gene2 = genes[2]
    
    
    data.frame(
      Genotype = rep(names(correlations), each =1),
      Correlations = sapply(names(correlations),function(genotype){
        correlations[[genotype]]$correlations[gene1,gene2]
      }),
      Gene_Pair = paste(gene1,gene2,sep = '_')
    )
  })
  
  data_combined = do.call(rbind,data_list)
  
  p= ggplot(data_combined, aes(x=Genotype, y=Correlations, group=Gene_Pair, color=Gene_Pair))+
    geom_line() + geom_point()+
    ggtitle('Gene Pairs Correlation Across Genotype')+
    ylab('Correlation')+
    theme_minimal()
  
  return(p)
  
}

plot_selected_pair = function(selected_pairs, correlations) {
  # Create a mapping for gene pair names to more readable names (e.g., ENSG IDs to gene names)
  gene_pair_names = c(
    "ENSG00000087258_ENSG00000130203" = "APOE-GNAO1",
    "ENSG00000130203_ENSG00000176340" = "APOE-COX8A",  # Replace "OtherGene1" with the actual gene name
    # Add more mappings as needed for other pairs, e.g., for APOE_neg_pair or additional pairs
    "ENSG00000111237_ENSG00000130203" = "APOE-VPS29",
    "ENSG00000130203_ENSG00000139597" = "APOE-N4BP2L1",
    "ENSG00000130203_ENSG00000166471" = "APOE-TMEM41B"
  )
  
  # Prepare the data for plotting
  data_list = lapply(selected_pairs, function(gene_pair) {
    genes = unlist(strsplit(gene_pair, '_'))
    gene1 = genes[1]
    gene2 = genes[2]
    
    data.frame(
      Genotype = rep(names(correlations), each = 1),
      Correlations = sapply(names(correlations), function(genotype) {
        correlations[[genotype]]$correlations[gene1, gene2]
      }),
      Gene_Pair = gene_pair  # Keep the original pair for mapping
    )
  })
  
  data_combined = do.call(rbind, data_list)
  
  # Map the gene pair names to readable names for the legend
  data_combined$Gene_Pair_Label = gene_pair_names[data_combined$Gene_Pair]
  
  # Create the plot with a modern, clean theme
  library(ggplot2)  # Ensure ggplot2 is loaded
  library(ggthemes)  # For additional themes (install if needed)
  library(ggrepel)
  
  p = ggplot(data_combined, aes(x = Genotype, y = Correlations, group = Gene_Pair_Label, color = Gene_Pair_Label)) +
    geom_line(linewidth = 1.5, alpha = 0.7) +  # Thicker, semi-transparent lines
    geom_point(size = 3) +  # Larger points
    geom_text_repel(aes(label = sprintf("%.2f", Correlations)),  # Use ggrepel for non-overlapping labels
                    size = 4,  # Increase text size for better visibility
                    color = "black",  # Ensure text is black for contrast
                    box.padding = 0.5,  # Add padding around labels
                    point.padding = 0.5) +  # Add padding around points
    ggtitle('Gene Pairs Correlation Across Genotype') +
    ylab('Correlation') +
    xlab('Genotype') +  # Add x-axis label
    theme_clean() +  # Use a modern, minimal theme
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Centered, bold title
      axis.title = element_text(size = 12),  # Larger axis titles
      axis.text = element_text(size = 10),  # Larger axis text
      legend.title = element_blank(),  # Remove legend title
      legend.position = "top",  # Move legend to top
      legend.text = element_text(size = 10),  # Larger legend text
      legend.background = element_blank(),  # Remove legend background
      legend.box.background = element_blank(),  # Remove legend box
      panel.grid.major = element_line(color = "gray95", linewidth = 0.3),  # Very light grid lines
      panel.grid.minor = element_blank()  # Remove minor grid lines
    ) +
    scale_color_brewer(palette = "Set2")  # Modern color palette
  
  return(p)
}


plot_selected_pair(APOE_pos_pair,correlations)
plot_selected_pair(APOE_neg_pair,correlations)

plot_selected_pair(COX8A_pos_pair,correlations)
plot_selected_pair(COX8A_neg_pair,correlations)




#Grafik ??izme
library(ggplot2)

plot_top_pairs = function(top_pairs, correlations_results) {
  plots = list()
  
  for(pattern_name in names(top_pairs)){
    pairs = top_pairs[[pattern_name]]
    
    if(nrow(pairs)>0 ){ #i??inde eleman var m??
      data_list = lapply(1:nrow(pairs), function(i){
        gene_pair = pairs$gene_pair[i]
        genes = unlist(strsplit(gene_pair, '_'))
        gene1 = genes[1]
        gene2 = genes[2]
        
        data.frame(
          Genotype = rep(names(correlations_results), each=1), #genotip isimlerini al??yorz
          Correlation = sapply(names(correlations_results), function(genotype){
            correlations_results[[genotype]]$correlations[gene1,gene2] # o genlerin cor de??erleri
          }),
          Gene_Pair = paste(gene1, gene2, sep = '_')
        )
        
      })
      
      data_combined = do.call(rbind, data_list)
      
      p = ggplot(data_combined, aes(x=Genotype, y=Correlation, group=Gene_Pair, color=Gene_Pair)) +
        geom_line()+ geom_point() +
        ggtitle(paste('Pattern:', pattern_name))+
        ylab('Correlation')+
        theme_minimal()
      
      plots[[pattern_name]] = p
    }
  }
  
  return(plots)
}

plots = plot_top_pairs(top_pairs, correlations)

#plotlar

for(plot_name in names(plots)){
  print(plots[plot_name])
}

count_unique_genes = function(filtered_patterns){
  unique_genes_count = list()
  
  for(pattern_name in names(filtered_patterns)){
    
    pattern = filtered_patterns[[pattern_name]]
    
    if(length(pattern)>0){
      gene_pairs = names(pattern)
      
      genes = unique(unlist(strsplit(gene_pairs,'_')))
      
      unique_genes_count[[pattern_name]] = list(count = length(genes),
                                                genes = genes)
    }
    
  }
  
  return(unique_genes_count)
}

unique_genes_count = count_unique_genes(filtered0.95)
unique_genes_count0.35 = count_unique_genes(filtered_patterns0.35)

library(biomaRt)
ensembl = useMart('ensembl', dataset = 'hsapiens_gene_ensembl')

merged_unique_symbol = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
            filters = 'ensembl_gene_id',
            values = unique_merged_genes_count,
            mart = ensembl)

ortak_symbol = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                           filters = 'ensembl_gene_id',
                           values = ortak,
                           mart = ensembl)








#largest clusterlardaki gen ??iftleri alma

neg_decreasing_largestCluster = read.csv('neg_cor_deacreasing_LargestCluster.csv')
neg_decreasing_largestCluster = neg_decreasing_largestCluster[,2]

pos_decreasing_largestCluster2 = read.csv('pos_cor_deacreasing_LargestCluster2.csv') 
pos_decreasing_largestCluster2 = pos_decreasing_largestCluster2[,2]

pos_decreasing_largestCluster1= read.csv('pos_cor_deacreasing_LargestCluster1.csv')
pos_decreasing_largestCluster1 = pos_decreasing_largestCluster1[,2]


#negdeki ??iftler  
neg_decreasing_largestCluster_genePairs = c() 
for(gene_pair in names(filteredPattern[['negative_cor_decreasing']])){
  
  genes_in_pair = unlist(strsplit(gene_pair, '_'))
  
  if(any(genes_in_pair %in% neg_decreasing_largestCluster)){
    neg_decreasing_largestCluster_genePairs = c(neg_decreasing_largestCluster_genePairs, gene_pair)
  }
  
}
length(neg_decreasing_largestCluster_genePairs) # cytoscape de kontrol edildi 4031 .



#pos 2 deki ??iftler
pos_decreasing_largestCluster2_genePairs = c()
for(gene_pair in names(filteredPattern[['positive_cor_decreasing']])){
  
  genes_in_pair = unlist(strsplit(gene_pair, '_'))
  
  if(any(genes_in_pair %in% pos_decreasing_largestCluster2)){
    pos_decreasing_largestCluster2_genePairs = c(pos_decreasing_largestCluster2_genePairs, gene_pair)
  }
  
}
length(pos_decreasing_largestCluster2_genePairs) # cytoscape de kontrol edildi 2. cluster 1832.


#pos 1 deki ??iftler  #2850 ??ift

pos_decreasing_largestCluster1_genePairs = c()
for(gene_pair in names(filteredPattern[['positive_cor_decreasing']])){
  
  genes_in_pair = unlist(strsplit(gene_pair, '_'))
  
  if(any(genes_in_pair %in% pos_decreasing_largestCluster1)){
    pos_decreasing_largestCluster1_genePairs = c(pos_decreasing_largestCluster1_genePairs, gene_pair)
  }
  
}




gene_pairs_pos = strsplit(pos_decreasing_largestCluster2_genePairs, '_')
df_pos_largestCluster2 = do.call(rbind, lapply(gene_pairs_pos, function(x) data.frame(Gene1=x[1],Gene2=x[2])))

#write.csv(df_pos_largestCluster2, 'df_pos_largestCluster2_genePairs.csv', row.names = FALSE)

gene_pairs_neg = strsplit(neg_decreasing_largestCluster_genePairs, '_')
df_neg_largestCluster= do.call(rbind, lapply(gene_pairs_neg, function(x) data.frame(Gene1=x[1],Gene2=x[2])))

#write.csv(df_neg_largestCluster, 'df_neg_largestCluster_genePairs.csv', row.names = F)


gene_pairs_posCluster1 = strsplit(pos_decreasing_largestCluster1_genePairs, '_')
df_pos_largestCluster1 = do.call(rbind, lapply(gene_pairs_posCluster1, function(x) data.frame(Gene1=x[1],Gene2=x[2])))

#write.csv(df_pos_largestCluster1, 'df_pos_largestCluster1_genePairs.csv', row.names = F)


# first ve second connections - APOE
# negative network
neg_cluster_gene_pairs = read.csv('df_neg_largestCluster_genePairs.csv')
pos_cluster2_gene_pairs = read.csv('df_pos_largestCluster2_genePairs.csv')


# x derece connection bulma
find_all_connections = function(df, gene, degree_limit =Inf){
  
  #connection listesi
  all_connections = data.frame(Gene1 = character(), Gene2 =character(), stringsAsFactors = F)
  #istenilen gen ile baslama
  current_genes = c(gene)
  degree = 1
  seen_genes = c(gene)
  
  while (length(current_genes) > 0 & degree <= degree_limit ) {
    #current gene i??eren t??m pairleri bulma
    new_connections = df[df$Gene1 %in% current_genes | df$Gene2 %in% current_genes,]
    
    #t??m connectionlara yenileri ekleme
    all_connections = unique(rbind(all_connections,new_connections))
    
    #bu patchteki connectionlar?? tespit etme
    next_genes = unique(c(new_connections$Gene1, new_connections$Gene2))
    
    #g??r??lm???? genleri sil
    current_genes = setdiff(next_genes, seen_genes)
    
    #yeni kesfedilen genleri ekle
    seen_genes = unique(c(seen_genes, current_genes))
    
    degree = degree + 1 
  }
  return(all_connections)
}

apoE_first_second_connections = find_all_connections(neg_cluster_gene_pairs, 'ENSG00000130203', 2) #neg networkteki
apoE_first_second_connections_pos2 = find_all_connections(pos_cluster2_gene_pairs, 'ENSG00000130203', 2) #pos networkteki




