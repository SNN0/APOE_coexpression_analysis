
setwd('N:/sinan/data2')

msbb_data=readRDS('MSBB_data.rds')

counts = msbb_data$counts
meta.data = msbb_data$meta_data

names(meta.data)[9] = 'apoe_genotype'

apoe_categ = c("23" = "E2E3", "33" = "E3E3", "34" = "E3E4", "22" = "E2E2", "24" = "E2E4", "44" = "E4E4")

meta.data = meta.data %>% mutate(
  apoe_genotype = factor(apoe_genotype, levels=names(apoe_categ), labels = apoe_categ)
)


dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = meta.data,
                             design = ~1)

dds = DESeq(dds)

norm.count = counts(dds, normalized = T)


dds.vst = vst(dds) #log2 
norm.count.vst = assay(dds.vst) %>% t()

long_non_coding=read.delim('N:/sinan/RNA_long_non-coding.txt', header = T, sep = '\t',stringsAsFactors = F)
long_non_coding=long_non_coding$ensembl_gene_id

coding_gene=read.delim('N:/sinan/gene_with_protein_product.txt', header = T, sep = '\t',stringsAsFactors = F)
coding_gene = coding_gene$ensembl_gene_id

mssb_gene_names = rownames(counts)
mssb_gene_names = sub("\\..*$", "", mssb_gene_names)

length(intersect(mssb_gene_names,coding_gene))
length(intersect(mssb_gene_names,long_non_coding))



library(foreach)
library(doParallel)

n_cores = 4  #max 32
cl = makeCluster(n_cores)
registerDoParallel(cl)

genotypes = c('E2E3', 'E3E3', 'E3E4')
correlations = cal_corr(norm.count.vst, genotypes)

stopCluster(cl)

library(CimpleG)

#save_object(correlations, 'MSBB_Correlations_padj.rds')

#load_object('MSBB_Correlations_padj.rds')

result2$common_pairs %>% filter(gene1 =='ENSG00000139597' |gene2=='ENSG00000139597')

gene_mapping <- c(
  "ENSG00000130203" = "APOE",
  "ENSG00000176340" = "COX8A",
  "ENSG00000087258" = "GNAO1",
  "ENSG00000111237" = "VPS29",
  "ENSG00000166471" = "TMEM41B",
  "ENSG00000139597" = "N4BP2L1"
)

row_n='ENSG00000111237'
col_n='ENSG00000029363'

row_m = grep(paste0('^',row_n,'\\.'), rownames(correlations$E2E3$correlations), value = T)
col_m = grep(paste0('^',col_n,'\\.'), colnames(correlations$E2E3$correlations), value = T)


correlations$E3E4$correlations[row_m, col_m]
correlations$E2E3$p_value[row_m, col_m]

#APOE-GNOA1 exp slope azaliyor. pval yok. cor azaliyor kriter yok
#APOE-N4B2PL1  exp slope azaliyor. pval yok.cor azaliyor kriter yok. y??ksek azal??yor.



#export from corr_sinan.R
plot_gene_correlation(meta.data, as.data.frame(norm.count.vst), "ENSG00000111237.19", "ENSG00000029363.17")






signi_pairs_msbb_e2e3 = filter_signi_corr_optimized(correlations, 'E2E3', threshold_corr = 0.70, threshold_p = 1)
pattern_msbb_0.05_dec = filter_correlation_pattern_all_optimized(signi_pairs_msbb_e2e3, correlations, threshold_diff = 0.1)


signi_pairs_msbb_e3e4 = filter_signi_corr_optimized(correlations, 'E3E4', threshold_corr = 0.70, threshold_p = 1)
pattern_msbb_0.05_inc = filter_correlation_pattern_increasing_E3E4(signi_pairs_msbb_e3e4, correlations, threshold_diff = 0.1)


saveRDS(pattern_msbb_0.05_dec, 'pattern/msbb_0.90cor_decreasing.rds')
saveRDS(pattern_msbb_0.05_dec, 'pattern/msbb_0.95cor_decreasing.rds')

saveRDS(pattern_msbb_0.05_inc, 'pattern/msbb_0.80cor_increasing.rds')#msbb de e3e4 te 0.05 den k??c??k 2 tane var.





pos_dec_MSBB = pattern_msbb_0.05_dec$positive_cor_decreasing
pos_dec_MSBB$gene1 <- sub("\\..*", "", pos_dec_MSBB$gene1)
pos_dec_MSBB$gene2 <- sub("\\..*", "", pos_dec_MSBB$gene2)

neg_dec_MSBB =pattern_msbb_0.05_dec$negative_cor_decreasing
neg_dec_MSBB$gene1 <- sub("\\..*", "", neg_dec_MSBB$gene1)
neg_dec_MSBB$gene2 <- sub("\\..*", "", neg_dec_MSBB$gene2)
msbb_dec_network = rbind(pos_dec_MSBB,neg_dec_MSBB)


pos_inc_MSBB = pattern_msbb_0.05_inc$positive_cor_increasing
pos_inc_MSBB$gene1 <- sub("\\..*", "", pos_inc_MSBB$gene1)
pos_inc_MSBB$gene2 <- sub("\\..*", "", pos_inc_MSBB$gene2)

neg_inc_MSBB =pattern_msbb_0.05_inc$negative_cor_increasing
neg_inc_MSBB$gene1 <- sub("\\..*", "", neg_inc_MSBB$gene1)
neg_inc_MSBB$gene2 <- sub("\\..*", "", neg_inc_MSBB$gene2)
msbb_inc_network = rbind(pos_inc_MSBB,neg_inc_MSBB)


rosmap_dec_network = read.csv('N:/sinan/pattern/rosmap_095_Decreasing_Network.csv')
rosmap_dec_network=as.data.table(rosmap_dec_network)

rosmap_inc_pattern = readRDS('N:/sinan/pattern/Rosmap_0.80cor_increasing.rds')
rosmap_pos_inc = rosmap_inc_pattern$positive_cor_increasing
rosmap_neg_inc = rosmap_inc_pattern$negative_cor_increasing
rosmap_inc_network = rbind(rosmap_pos_inc,rosmap_neg_inc)


second_degree_APOE_network = readxl::read_xlsx('N:/sinan/second_con_APOE.xlsx')
second_degree_APOE_network = as.data.table(second_degree_APOE_network)


find_common_pairs = function(data1, data2) {
  # Gen C'iftlerinin sD1rasD1nD1 dC<zenleyerek her iki veri setinde de aynD1 olacak Eekilde ayarlD1yoruz
  reorder_pairs = function(data) {
    # Her bir satD1rda gen C'iftlerini alfabetik sD1raya gC6re dizdim
    data[, .(gene1 = pmin(gene1, gene2), gene2 = pmax(gene1, gene2))]
  }
  
  # D0ki veri setindeki gen C'iftlerini sD1ralayD1p dC<zenliyoruz
  ordered_data1 = reorder_pairs(data1)
  ordered_data2 = reorder_pairs(data2)
  
  # Ortak gen C'iftlerini bulmak iC'in birleEtiriyoruz (inner join)
  common_pairs = merge(ordered_data1, ordered_data2, by = c("gene1", "gene2"))
  
  
  list(common_pairs = common_pairs, count_common = nrow(common_pairs))
}


result3 = find_common_pairs(msbb_dec_network, second_degree_APOE_network)


result2 = find_common_pairs(rosmap_dec_network, msbb_dec_network)
result2$common_pairs %>% filter(gene1 =='ENSG00000139597' |gene2=='ENSG00000139597')



result = find_common_pairs(rosmap_inc_network, msbb_inc_network)
result$common_pairs %>% filter(gene1 =='ENSG00000130203' |gene2=='ENSG00000130203')




#ORTAK --  min cor 0.7 , dif 0.1, pval =1
#COX8A - DERL1 - sari node, pval 0.07
#COX8A - KPNA3
#GNOA1 - RBM27
#GNOA1 - JAK2
#GNOA1 - RNF180
#GNOA1 - SHARPIN
#GNOA1 - TRIQK
#VPS29 - BCLAF1 - sari node  pval 0.06
#TMEM41B - RPAP3
#TMEM41B - BCLAF1 - sari node  pval 0.19
#TMEM41B - AP5M1
#N4BP2L1 - CST3 - sari node  pval 0.13
#N4BP2L1 - CAMSAP2
#N4BP2L1 - ORMDL1 


#----------------------MSBB-Transcriptome vs ROSMAP-Transcriptome---------------------------------
#ORTAK --  min cor 0.7 , dif 0.1, raw-pval =0.05
#COX8A -KPNA3
#COX8A - DERL1 - sari node
#GNOA1 - RBM27
#GNOA1 - JAK2
#GNOA1 - RNF180
#GNOA1 - SHARPIN
#GNOA1 - TRIQK
#VPS29 - BCLAF1 - sari node 
#TMEM41B - RPAP3
#TMEM41B - BCLAF1 - sari node 
#TMEM41B - AP5M1
#N4BP2L1 - CST3 - sari node  
#N4BP2L1 - CAMSAP2
#N4BP2L1 - ORMDL1 
#-------------------------------------------------------



result2 = find_common_pairs(neg_dec_MSBB, neg_dec_rosmap)


pos_dec_MSBB.unique = unique(c(pos_dec_MSBB$gene1,pos_dec_MSBB$gene2))
pos_dec_rosmap.unique = unique(c(pos_dec_rosmap$gene1,pos_dec_rosmap$gene2))

neg_dec_MSBB.unique = unique(c(neg_dec_MSBB$gene1,neg_dec_MSBB$gene2))
neg_dec_rosmap.unique = unique(c(neg_dec_rosmap$gene1,neg_dec_rosmap$gene2))



length(intersect(pos_dec_MSBB.unique,pos_dec_rosmap.unique)) #1109 ortak gen

length(intersect(neg_dec_MSBB.unique,neg_dec_rosmap.unique)) #766 ortak gen
