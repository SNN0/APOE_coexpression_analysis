
#----------------------------------------------




setwd('C:/Users/sinan/Desktop/staj_calisma/patternYeni')

rosmap_transcriptome_decreasing = readRDS('Rosmap_0.95cor_decreasing.rds')
rosmap_transcriptome_increasing = readRDS('Rosmap_0.80cor_increasing.rds')

#writexl::write_xlsx(rosmap_transcriptome_decreasing$negative_cor_decreasing, 'rosmap_090cor_decreasing_neg.xlsx')


msbb_transcriptome_decreasing = readRDS('msbb_0.95cor_decreasing.rds')


rosmap_proteom_decreasing = readRDS('rosmap_proteom_0.95cor_decreasing.rds')
rosmap_proteom_increasing = readRDS('rosmap_proteom_0.80cor_increasing.rds')


#RosMap Transcriptome pattern
rosmap_tc_pos_decreasing = rosmap_transcriptome_decreasing$positive_cor_decreasing
rosmap_tc_neg_decreasing = rosmap_transcriptome_decreasing$negative_cor_decreasing
rosmap_tc_pos_increasing = rosmap_transcriptome_increasing$positive_cor_increasing
rosmap_tc_neg_increasing = rosmap_transcriptome_increasing$negative_cor_increasing


#MSBB Transcriptoe pattern
msbb_tc_pos_decreasing = msbb_transcriptome_decreasing$positive_cor_decreasing
msbb_tc_neg_decreasing = msbb_transcriptome_decreasing$negative_cor_decreasing
#msbb_tc_pos_increasing = msbb_transcriptome_increasing$positive_cor_increasing
#msbb_tc_neg_increasing = msbb_transcriptome_increasing$negative_cor_increasing

msbb_tc_pos_decreasing$gene1 <- sub("\\..*", "", msbb_tc_pos_decreasing$gene1)
msbb_tc_pos_decreasing$gene2 <- sub("\\..*", "", msbb_tc_pos_decreasing$gene2)

msbb_tc_neg_decreasing$gene1 <- sub("\\..*", "", msbb_tc_neg_decreasing$gene1)
msbb_tc_neg_decreasing$gene2 <- sub("\\..*", "", msbb_tc_neg_decreasing$gene2)

#msbb_tc_pos_increasing$gene1 <- sub("\\..*", "", msbb_tc_pos_increasing$gene1)
#msbb_tc_pos_increasing$gene2 <- sub("\\..*", "", msbb_tc_pos_increasing$gene2)

#msbb_tc_neg_increasing$gene1 <- sub("\\..*", "", msbb_tc_neg_increasing$gene1)
#msbb_tc_neg_increasing$gene2 <- sub("\\..*", "", msbb_tc_neg_increasing$gene2)


#ROSMAP Proteom pattern

rosmap_proteom_pos_decreasing = rosmap_proteom_decreasing$positive_cor_decreasing
rosmap_proteom_neg_decreasing = rosmap_proteom_decreasing$negative_cor_decreasing
rosmap_proteom_pos_increasing = rosmap_proteom_increasing$positive_cor_increasing
rosmap_proteom_neg_increasing = rosmap_proteom_increasing$negative_cor_increasing

rosmap_proteom_pos_decreasing$gene1 <- sub("\\..*", "", rosmap_proteom_pos_decreasing$gene1)
rosmap_proteom_pos_decreasing$gene2 <- sub("\\..*", "", rosmap_proteom_pos_decreasing$gene2)

rosmap_proteom_neg_decreasing$gene1 <- sub("\\..*", "", rosmap_proteom_neg_decreasing$gene1)
rosmap_proteom_neg_decreasing$gene2 <- sub("\\..*", "", rosmap_proteom_neg_decreasing$gene2)

rosmap_proteom_pos_increasing$gene1 <- sub("\\..*", "", rosmap_proteom_pos_increasing$gene1)
rosmap_proteom_pos_increasing$gene2 <- sub("\\..*", "", rosmap_proteom_pos_increasing$gene2)

rosmap_proteom_neg_increasing$gene1 <- sub("\\..*", "", rosmap_proteom_neg_increasing$gene1)
rosmap_proteom_neg_increasing$gene2 <- sub("\\..*", "", rosmap_proteom_neg_increasing$gene2)

# ENSMBL donusturme
rosmap_proteom_pos_decreasing = convert_symbol_to_ensembl(rosmap_proteom_pos_decreasing, 'gene1','gene2')
rosmap_proteom_neg_decreasing = convert_symbol_to_ensembl(rosmap_proteom_neg_decreasing, 'gene1','gene2')
rosmap_proteom_pos_increasing = convert_symbol_to_ensembl(rosmap_proteom_pos_increasing, 'gene1','gene2')
rosmap_proteom_neg_increasing = convert_symbol_to_ensembl(rosmap_proteom_neg_increasing, 'gene1','gene2')

rosmap_proteom_pos_decreasing$gene1 = rosmap_proteom_pos_decreasing$ensembl_gene_id_gene1
rosmap_proteom_pos_decreasing$gene2 = rosmap_proteom_pos_decreasing$ensembl_gene_id_gene2

rosmap_proteom_neg_decreasing$gene1 = rosmap_proteom_neg_decreasing$ensembl_gene_id_gene1
rosmap_proteom_neg_decreasing$gene2 = rosmap_proteom_neg_decreasing$ensembl_gene_id_gene2

rosmap_proteom_pos_increasing$gene1 = rosmap_proteom_pos_increasing$ensembl_gene_id_gene1
rosmap_proteom_pos_increasing$gene2 = rosmap_proteom_pos_increasing$ensembl_gene_id_gene2

rosmap_proteom_neg_increasing$gene1 = rosmap_proteom_neg_increasing$ensembl_gene_id_gene1
rosmap_proteom_neg_increasing$gene2 = rosmap_proteom_neg_increasing$ensembl_gene_id_gene2


# func --*******************----------------------------------

find_common_pairs = function(data1, data2) {
        # Gen çiftlerinin sırasını düzenleyerek her iki veri setinde de aynı olacak şekilde ayarlıyoruz
        reorder_pairs = function(data) {
                # Her bir satırda gen çiftlerini alfabetik sıraya göre dizdim
                data[, .(gene1 = pmin(gene1, gene2), gene2 = pmax(gene1, gene2))]
        }
        
        # İki veri setindeki gen çiftlerini sıralayıp düzenliyoruz
        ordered_data1 = reorder_pairs(data1)
        ordered_data2 = reorder_pairs(data2)
        
        # Ortak gen çiftlerini bulmak için birleştiriyoruz (inner join)
        common_pairs = merge(ordered_data1, ordered_data2, by = c("gene1", "gene2"))
        
        
        list(common_pairs = common_pairs, count_common = nrow(common_pairs))
}

library(biomaRt)
convert_symbol_to_ensembl <- function(df, gene1_col, gene2_col) {
        
        # Ensembl mart'ı kullanarak insan genomunu seçelim
        ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        
        # Gene symbol'leri toplayalım
        gene_symbols <- unique(c(df[[gene1_col]], df[[gene2_col]]))
        
        # Gene symbol'leri Ensembl ID'leriyle eşleştirelim
        mapping <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id'),
                         filters = 'hgnc_symbol', 
                         values = gene_symbols, 
                         mart = ensembl)
        
        # Boş (bulunmayan) sembolleri NA yapalım
        mapping$ensembl_gene_id[mapping$ensembl_gene_id == ""] <- NA
        
        # İlk gen kolonu için eşleşmeleri yapalım
        df <- merge(df, mapping, by.x = gene1_col, by.y = 'hgnc_symbol', all.x = TRUE)
        
        # İkinci gen kolonu için eşleşmeleri yapalım
        df <- merge(df, mapping, by.x = gene2_col, by.y = 'hgnc_symbol', all.x = TRUE, suffixes = c('_gene1', '_gene2'))
        
        # NA olan değerleri yerinde bırakalım
        df$ensembl_gene_id_gene1[is.na(df$ensembl_gene_id_gene1)] <- NA
        df$ensembl_gene_id_gene2[is.na(df$ensembl_gene_id_gene2)] <- NA
        
        return(df)
}


# kaç ortak var A-B  ikilisi var  ?? ?! ----------------------

rosmap_transcriptome_vs_rosmap_proteom = find_common_pairs(rosmap_proteom_pos_decreasing, rosmap_tc_pos_decreasing)
rosmap_transcriptome_vs_rosmap_proteom = find_common_pairs(rosmap_proteom_neg_decreasing, rosmap_tc_neg_decreasing)
rosmap_transcriptome_vs_rosmap_proteom = find_common_pairs(rosmap_proteom_pos_increasing, rosmap_tc_pos_increasing)
rosmap_transcriptome_vs_rosmap_proteom = find_common_pairs(rosmap_proteom_neg_increasing, rosmap_tc_neg_increasing)

rosmap_transcriptome_vs_msbb_transcriptome = find_common_pairs(rosmap_tc_neg_decreasing, msbb_tc_neg_decreasing)
rosmap_transcriptome_vs_msbb_transcriptome = find_common_pairs(rosmap_tc_pos_decreasing, msbb_tc_pos_decreasing)
rosmap_transcriptome_vs_msbb_transcriptome = find_common_pairs(rosmap_tc_neg_increasing, msbb_tc_neg_increasing)
rosmap_transcriptome_vs_msbb_transcriptome = find_common_pairs(rosmap_tc_pos_increasing, msbb_tc_pos_increasing)

msbb_transcriptome_vs_rosmap_proteom = find_common_pairs(msbb_tc_neg_decreasing,rosmap_proteom_neg_decreasing)
msbb_transcriptome_vs_rosmap_proteom = find_common_pairs(msbb_tc_pos_decreasing,rosmap_proteom_pos_decreasing)
msbb_transcriptome_vs_rosmap_proteom = find_common_pairs(msbb_tc_neg_increasing,rosmap_proteom_neg_increasing)
msbb_transcriptome_vs_rosmap_proteom = find_common_pairs(msbb_tc_pos_increasing,rosmap_proteom_pos_increasing)


# rosmap_trans vs msbb_trans - decreasing pos ta  ortak 5 tane var. !!!!!!!!!!!!!!!


#-*--------------------------------------------------------------**----------------

# direkt ortak olarak kaç gen,protein var?.  Unique gen her pattern için ve kaç ortak

# unique sayıları
length(unique(c(rosmap_proteom_pos_increasing$gene1, rosmap_proteom_pos_increasing$gene2)))

#unique Listeler

rosmap_tc_neg_decreasing_unique = unique(c(rosmap_tc_neg_decreasing$gene1, rosmap_tc_neg_decreasing$gene2))
rosmap_tc_neg_increasing_unique = unique(c(rosmap_tc_neg_increasing$gene1, rosmap_tc_neg_increasing$gene2))
rosmap_tc_pos_increasing_unique = unique(c(rosmap_tc_pos_increasing$gene1, rosmap_tc_pos_increasing$gene2))
rosmap_tc_pos_decreasing_unique = unique(c(rosmap_tc_pos_decreasing$gene1, rosmap_tc_pos_decreasing$gene2))

msbb_tc_neg_decreasing_unique = unique(c(msbb_tc_neg_decreasing$gene1, msbb_tc_neg_decreasing$gene2))
msbb_tc_pos_decreasing_unique = unique(c(msbb_tc_pos_decreasing$gene1, msbb_tc_pos_decreasing$gene2))
msbb_tc_pos_increasing_unique = unique(c(msbb_tc_pos_increasing$gene1, msbb_tc_pos_increasing$gene2))
msbb_tc_neg_increasing_unique = unique(c(msbb_tc_neg_increasing$gene1, msbb_tc_neg_increasing$gene2))

rosmap_proteom_neg_decreasing_unique = unique(c(rosmap_proteom_neg_decreasing$gene1, rosmap_proteom_neg_decreasing$gene2))
rosmap_proteom_pos_decreasing_unique = unique(c(rosmap_proteom_pos_decreasing$gene1, rosmap_proteom_pos_decreasing$gene2))
rosmap_proteom_pos_increasing_unique = unique(c(rosmap_proteom_pos_increasing$gene1, rosmap_proteom_pos_increasing$gene2))
rosmap_proteom_neg_increasing_unique = unique(c(rosmap_proteom_neg_increasing$gene1, rosmap_proteom_neg_increasing$gene2))


#----------------------------
# ENSG00000130203 - APOE

#datasetler arasında ortak sayı

#rosmap-trans  vs rosmap-proteom
rosmap_transcriptome_vs_rosmap_proteom_pos_dec= intersect(rosmap_tc_pos_decreasing_unique, rosmap_proteom_pos_decreasing_unique)
intersect(rosmap_transcriptome_vs_rosmap_proteom_pos_dec, msbb_tc_pos_decreasing_unique)


# inc ve dec leri birleştirme.  pos neg diye bakmıyorum ikiside aynı
rosmap_trans_dec = rbind(rosmap_tc_neg_decreasing, rosmap_tc_pos_decreasing)
rosmap_trans_inc = rbind(rosmap_tc_neg_increasing, rosmap_tc_pos_increasing)

msbb_trans_dec = rbind(msbb_tc_neg_decreasing,msbb_tc_pos_decreasing)
msbb_trans_inc = rbind(msbb_tc_neg_increasing,msbb_tc_pos_increasing)

rosmap_proteom_dec = rbind(rosmap_proteom_neg_decreasing, rosmap_proteom_pos_decreasing)
rosmap_proteom_inc =  rbind(rosmap_proteom_neg_increasing, rosmap_proteom_pos_increasing)

#buradaki rosmap için olan largest network :
rosmap_trans_dec_network = read.csv('rosmap095Decreasing_NETWORK.csv')

rosmap_trans_dec_network$gene1 <- sub(" \\(.*", "", rosmap_trans_dec_network$name)
rosmap_trans_dec_network$gene2 <- sub(".*\\) ", "", rosmap_trans_dec_network$name)

library(data.table)
rosmap_trans_dec_network <- data.table(gene1 = rosmap_trans_dec_network$gene1, gene2 = rosmap_trans_dec_network$gene2)
write.csv(rosmap_trans_dec_network, 'rosmap_095_Decreasing_Network.csv' )


#buda msbb için olan largerst network:
msbb_trans_dec_network = read.csv('C:/Users/sinan/Documents/m.csv')

msbb_trans_dec_network$gene1 <- sub(" \\(.*", "", msbb_trans_dec_network$name)
msbb_trans_dec_network$gene2 <- sub(".*\\) ", "", msbb_trans_dec_network$name)

msbb_trans_dec_network <- data.table(gene1 = msbb_trans_dec_network$gene1, gene2 = msbb_trans_dec_network$gene2)

msbb_trans_dec_network$gene1 <- sub("\\..*", "", msbb_trans_dec_network$gene1)
msbb_trans_dec_network$gene2 <- sub("\\..*", "", msbb_trans_dec_network$gene2)


#----------
rosmap_proteom_rosmap_trans_msbb_trans_dec = rbind(rosmap_trans_dec_network, msbb_trans_dec_network)
rosmap_proteom_rosmap_trans_msbb_trans_dec = rbind(rosmap_proteom_rosmap_trans_msbb_trans_dec, rosmap_proteom_dec[,1:2])

rosmap_trans_rosmap_proteom_inc = rbind(rosmap_trans_inc,rosmap_proteom_inc[,1:3])


#--------------


rosmap_trans_dec_unique = unique(c(rosmap_trans_dec$gene1,rosmap_trans_dec$gene2)) #4989 - word de farklı çünkü bundan önce en büyük networkü seçtim.
rosmap_trans_inc_unique = unique(c(rosmap_trans_inc$gene1,rosmap_trans_inc$gene2))

#bu da networlü olan hali:
rosmap_trans_dec_network_unique =  unique(c(rosmap_trans_dec_network$gene1,rosmap_trans_dec_network$gene2)) #3394

msbb_trans_dec_unique = unique(c(msbb_trans_dec$gene1,msbb_trans_dec$gene2))
msbb_trans_inc_unique = unique(c(msbb_trans_inc$gene1,msbb_trans_inc$gene2))

#bu da networklü hali: 
msbb_trans_dec_network_unique = unique(c(msbb_trans_dec_network$gene1,msbb_trans_dec_network$gene2)) #2078


rosmap_proteom_dec_unique = unique(c(rosmap_proteom_dec$gene1,rosmap_proteom_dec$gene2))
rosmap_proteom_inc_unique = unique(c(rosmap_proteom_inc$gene1,rosmap_proteom_inc$gene2))

length(rosmap_proteom_inc_unique)

# birleşmiş (dec ve inc olarak) da ortak var mı  A-B şeklinde?

#rosmap_t vs msbb_t - dec  1!! 0 tane var. !!!  büyük networkten önce bu
find_common_pairs(rosmap_trans_dec, msbb_trans_dec) 

#rosmap_t vs msbb_t - dec  1!! büyük networkten sonra :0  tane var
find_common_pairs(rosmap_trans_dec_network, msbb_trans_dec_network) 

#rosmap-t vs msbb-t - inc  yok.
find_common_pairs(rosmap_trans_inc, msbb_trans_inc)

#rosmap-t vs rosmap-proteom dec yok.
find_common_pairs(rosmap_trans_dec, rosmap_proteom_dec)

#rosmap -t vs rosmap-proteom inc 
find_common_pairs(rosmap_trans_inc, rosmap_proteom_inc)

#rosmap-proteom vs msbb-t dec : yok
find_common_pairs(rosmap_proteom_dec, msbb_trans_dec)

#rosmap-proteom vs msbb-t inc : 

find_common_pairs(rosmap_proteom_inc, msbb_trans_inc)



#şimdi unique olarak ortaklara bakıcaz.

#rosmap-t vs msbb-t dec
length(intersect(rosmap_trans_dec_network_unique, msbb_trans_dec_network_unique))
rosmap_t_vs_msbb_t_dec = intersect(rosmap_trans_dec_network_unique, msbb_trans_dec_network_unique)

#rosmap-t vs msbb-t inc
length(intersect(rosmap_trans_inc_unique, msbb_trans_inc_unique))
rosmap_t_vs_msbb_t_inc = intersect(rosmap_trans_inc_unique, msbb_trans_inc_unique)

#rosmap-t vs rosmap-proteom  dec
length(intersect(rosmap_trans_dec_network_unique, rosmap_proteom_dec_unique))
rosmap_t_vs_rosmap_p_dec = intersect(rosmap_trans_dec_network_unique, rosmap_proteom_dec_unique)

#rosmap-t vs rosmap-proteom  inc #13 tane
length(intersect(rosmap_trans_inc_unique, rosmap_proteom_inc_unique))
rosmap_t_vs_rosmap_p_inc= intersect(rosmap_trans_inc_unique, rosmap_proteom_inc_unique)


#msbb-t vs rosmap-proteom  dec
length(intersect(msbb_trans_dec_network_unique, rosmap_proteom_dec_unique))
rosmap_p_vs_msbb_t_dec = intersect(msbb_trans_dec_unique, rosmap_proteom_dec_unique)


#msbb-t vs rosmap-proteom  inc
length(intersect(msbb_trans_inc_unique, rosmap_proteom_inc_unique))


#Rosmap-t vs msbb-t vs rosmap-proteom dec:
length(intersect(intersect(rosmap_trans_dec_network_unique,msbb_trans_dec_network_unique),rosmap_proteom_dec_unique))
rosmap_t_vs_rosmap_p_vs_msbb_t_dec = intersect(intersect(rosmap_trans_dec_network_unique,msbb_trans_dec_network_unique),rosmap_proteom_dec_unique)
write.csv(rosmap_t_vs_rosmap_p_vs_msbb_t_dec, '3lukesisim_dec.csv',row.names = F)


#Rosmap-t vs msbb-t vs rosmap-proteom inc:
length(intersect(intersect(rosmap_trans_inc_unique,msbb_trans_inc_unique),rosmap_proteom_inc_unique))


#Enrichment Analysis :  ---------------------------------

library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)


enrichment = function(x,y,title){
        plot=enrichGO(
                x,
                org.Hs.eg.db,
                keyType = 'ENSEMBL',
                ont = y,
                pvalueCutoff = 0.05,
                pAdjustMethod = 'BH',
                qvalueCutoff = 0.2,
                minGSSize = 10,
                maxGSSize = 500,
                readable = TRUE,
                pool = FALSE
        )
        return(plot)
}

enrichment_KEGG = function(x,title){
        entrezID = unlist(mget(x, envir=org.Hs.egENSEMBL2EG,
                               ifnotfound = NA))
        
        plot=enrichKEGG(entrezID, pvalueCutoff = 0.05)
        return(plot)
}



#enrichment analysis *------------------------------------------
#rosmap_t_vs_msbb_t_inc
#rosmap_t_vs_rosmap_p_vs_msbb_t_dec

#rosmap_proteom_dec_unique sadece proteom dec 

#rosmap_t_vs_msbb_t_dec 591 olan
#rosmap_t_vs_rosmap_p_dec  92 olan

#genel olarak üstekiler incelendi fakat birde rosmap_t vs msbb_t  dec kesişimi inceleme ? 

azalan_kesisim_3lü = read.csv('azalan_kesisim_3lü.csv')
azalan_kesisim_3lü$converted_alias

artan_kesisim_2li = read.csv('artan_kesisim_2li.csv')
artan_kesisim_2li$converted_alias



sonuc=enrichment(na.omit(rosmap_t_vs_msbb_t_dec) , 'BP', 'pos_cor_decreasing2_GO-BP')
sonuc1=enrichment(na.omit(rosmap_t_vs_msbb_t_dec) , 'CC', 'pos_cor_decreasing2_GO-CC')
sonuc2=enrichment(na.omit(rosmap_t_vs_msbb_t_dec)  , 'MF', 'pos_cor_decreasing2_GO-MF')
sonucKEGG =enrichment_KEGG(na.omit(rosmap_t_vs_msbb_t_dec) , 'pos_cor_decresing2_KEGG')
sonucKEGG = select(sonucKEGG, -c('category','subcategory'))

all_results <- rbind(
        sonuc@result %>% mutate(Category = "BP"),
        sonuc1@result %>% mutate(Category = "CC"),
        sonuc2@result %>% mutate(Category = 'MF'),
        sonucKEGG@result %>% mutate(Category = "KEGG")
)

all_results <- all_results %>%
        mutate(p.adjust = as.numeric(p.adjust),
               GeneRatio = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))))


selected_descriptions <- c('Alzheimer disease','Endocytosis','sphingomyelin biosynthetic process',
                           'vesicle organization',
                           'neuron to neuron synapse', 'asymmetric synapse',
                           'postsynaptic density', 'Ubiquitin mediated proteolysis','protein polyubiquitination',
                           'dendrite development','transport along microtubule','macroautophagy','aerobic respiration',
                           'GTPase activity'
)

selected_descriptions_2 <- c('Alzheimer disease','Oxidative phosphorylation','aerobic respiration',
                           'respiratory chain complex',
                           'purine nucleoside triphosphate biosynthetic process', 'ATP biosynthetic process',
                           'NADH dehydrogenase (ubiquinone) activity', 'proton motive force-driven mitochondrial ATP synthesis','mitochondrial protein-containing complex',
                           'mitochondrial ATP synthesis coupled electron transport'
)

selected_descriptions_3 <- c('Alzheimer disease','Amyotrophic lateral sclerosis','macroautophagy',
                           'vesicle organization',
                           'microtubule organizing center organization', 'sphingomyelin biosynthetic process',
                           'protein polyubiquitination', 'regulation of postsynapse organization','phospholipid biosynthetic process',
                           'calcium ion transport'
)

# Veriyi filtreleme
filtered_results <- all_results %>%
        filter(Description %in% selected_descriptions_3)


# Veriyi filtreleme
filtered_results_2 <- all_results %>%
        filter(Description %in% selected_descriptions_2)


# 1) Pattern değişkeni ekle ve birleştir
filtered_results$Pattern <- "Decreasing Pattern"
filtered_results_2$Pattern <- "Increasing Pattern"


combined_data <- rbind(filtered_results, filtered_results_2)

# -log10(p-value) sütunu
combined_data$neg_log10_pvalue <- -log10(combined_data$pvalue)

# GO:ID kısmını at
combined_data$Description <- gsub("GO:[0-9]+ ", "", combined_data$Description)

# Pattern’i factor olarak ayarla, facet sırası bozulmasın:
combined_data$Pattern <- factor(combined_data$Pattern,
                                levels = c("Decreasing Pattern","Increasing Pattern"))


my_order <- c(
        "Alzheimer disease",
        "aerobic respiration",
        "sphingomyelin biosynthetic process",
        "vesicle organization",
        "neuron to neuron synapse","asymmetric synapse",'postsynaptic density','Endocytosis',
        'Ubiquitin mediated proteolysis','protein polyubiquitination',
        'dendrite development','transport along microtubule','macroautophagy',
        'GTPase activity','Oxidative phosphorylation','respiratory chain complex',
        'purine nucleoside triphosphate biosynthetic process','ATP biosynthetic process',
        'NADH dehydrogenase (ubiquinone) activity','proton motive force-driven mitochondrial ATP synthesis',
        'mitochondrial protein-containing complex', 'mitochondrial ATP synthesis coupled electron transport'
        
)



combined_data$Description <- factor(
        combined_data$Description,
        levels = rev(my_order)
)

# --- 2) Plot ---
p <- ggplot(combined_data,
            aes(x = Pattern,
                y = Description)) +
        # Noktaları sadece size ile çiz, sabit kırmızı renk
        geom_point(aes(size = neg_log10_pvalue),
                   colour = "firebrick") +
        # Yan yana iki panel; x ölçeği özgür, y ortak
        facet_grid(. ~ Pattern,
                   scales = "free_x",
                   space  = "free_x",
                   labeller = labeller(
                           Pattern = c("Decreasing Pattern" = "Decreasing Pattern",
                                       "Increasing Pattern" = "Increasing Pattern")
                   )
        ) +
        # Sembol büyüklüğü legend’ı
        scale_size_continuous(name = "-log10(adj.pvalue)",
                              range = c(3, 8)) +
        # Eksensiz, başlıksız
        labs(x = NULL, y = NULL) +
        # Hafif grid’li, temiz arka plan
        theme_bw() +
        theme(
                # Facet başlıklarını hoş bir kutu içine al
                strip.background = element_rect(fill = "lightgrey", colour = "black"),
                strip.text       = element_text(face = "bold", size = 12),
                # Major grid çizgileri açık, minor kapalı
                panel.grid.major = element_line(color = "grey90"),
                panel.grid.minor = element_blank(),
                # Y-eksenindeki terimler biraz daha okunaklı
                axis.text.y      = element_text(size = 10),
                # X-eksenini tamamen gizle
                axis.text.x      = element_blank(),
                axis.ticks.x     = element_blank(),
                # İki panel arasındaki boşluk
                panel.spacing    = unit(1, "lines"),
                # Legend sağda
                legend.position  = "right",
                legend.title     = element_text(size = 10),
                legend.text      = element_text(size = 8)
        )

# Plot’u ekrana bas
print(p)


p <- ggplot(combined_data,
            aes(x = 1,
                y = reorder(Description, neg_log10_pvalue))) +
        geom_point(aes(size = neg_log10_pvalue),
                   colour = "firebrick") +
        scale_size_continuous(name = "-log10(adj.pvalue)",
                              range = c(3, 8)) +
        labs(x = NULL, y = NULL) +
        theme_bw() +
        theme(
                # Panel arka planını beyaz bırak
                panel.background    = element_rect(fill = "white"),
                # Hem yatay hem dikey major grid çizgileri
                panel.grid.major    = element_line(color = "grey80", size = 0.5),
                # Opsiyonel: daha ince minor grid eklemek istersen aç
                panel.grid.minor    = element_line(color = "grey90", size = 0.25),
                # X-ekseni yazı ve tıklardan kurtul
                axis.text.x         = element_blank(),
                axis.ticks.x        = element_blank(),
                # Y-etiketleri okunaklı
                axis.text.y         = element_text(size = 10),
                # Legend ayarları
                legend.position     = "right",
                legend.title        = element_text(size = 10),
                legend.text         = element_text(size = 8)
        ) +
        scale_x_continuous(expand = expansion(add = 0.2))



print(p)


library(ggplot2)
merged_enriched=ggplot(filtered_results, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
        geom_point() +
        scale_color_gradient(low = "red", high = "blue") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Gene Ratio", y = "Description", color = "p.adj", size = "Count")

ggsave(filename = 'dec_248_ad_related/enrichment_rosmap_msbb_proteom_dec_248.png' , plot = merged_enriched, device = "png")
ggsave(filename = 'azalan_kesisim_3lu_enrichment.svg' , plot = merged_enriched, device = "svg",width = 6.5,height = 4.5)





#AD related gene
AD_related_gene_list = read.csv('AD_related_gene_list.yeni.csv')

which(rosmap_t_vs_msbb_t_dec %in% AD_related_gene_list$ENSEMBL)
which(rosmap_t_vs_rosmap_p_dec %in% AD_related_gene_list$ENSEMBL)
which(rosmap_p_vs_msbb_t_dec %in% AD_related_gene_list$ENSEMBL)
which(rosmap_t_vs_rosmap_p_vs_msbb_t_dec %in% AD_related_gene_list$ENSEMBL)



which(rosmap_t_vs_msbb_t_dec %in% AD_related_gene_list$ENSEMBL) # 10 tane

rosmap_t_vs_msbb_t_adrelated=rosmap_t_vs_msbb_t_dec[which(rosmap_t_vs_msbb_t_dec %in% AD_related_gene_list$ENSEMBL)]
rosmap_t_vs_rosmap_p_adrelated=rosmap_t_vs_rosmap_p_dec[which(rosmap_t_vs_rosmap_p_dec %in% AD_related_gene_list$ENSEMBL)]
rosmap_p_vs_msbb_t_adrelated=rosmap_p_vs_msbb_t_dec[which(rosmap_p_vs_msbb_t_dec %in% AD_related_gene_list$ENSEMBL)]

ad_related_con_rosmapT_vsMSBBt_dec=find_all_connections_multiple_genes(rosmap_proteom_rosmap_trans_msbb_trans_dec, rosmap_t_vs_msbb_t_adrelated,2)
ad_related_con_rosmapT_vsROSMAPp_dec=find_all_connections_multiple_genes(rosmap_proteom_rosmap_trans_msbb_trans_dec, rosmap_t_vs_rosmap_p_adrelated,2)
ad_related_con_rosmapP_vsMSBBt_dec=find_all_connections_multiple_genes(rosmap_proteom_rosmap_trans_msbb_trans_dec, rosmap_p_vs_msbb_t_adrelated,2)


uclu_kesisim_dec_interaction = find_all_connections_multiple_genes(rosmap_proteom_rosmap_trans_msbb_trans_dec,rosmap_t_vs_rosmap_p_vs_msbb_t_dec,1)
write_all_dataframes_to_excel(uclu_kesisim_dec_interaction, 'uclu_kesisim_dec')


#rosmap_proteom_rosmap_trans_msbb_trans_dec[rosmap_proteom_rosmap_trans_msbb_trans_dec$gene2=='ENSG00000095564',]

write_all_dataframes_to_excel(ad_related_con_rosmapT_vsMSBBt_dec, 'romap_t_msbb_t_dec_v2')
write_all_dataframes_to_excel(ad_related_con_rosmapT_vsROSMAPp_dec, 'romap_t_rosmap_p_dec_v2')
write_all_dataframes_to_excel(ad_related_con_rosmapP_vsMSBBt_dec, 'romap_p_msbb_t_dec_v2')


#ad related gene  inc

rosmap_p_adrelated_inc=rosmap_proteom_inc_unique[which(rosmap_proteom_inc_unique %in% AD_related_gene_list$ENSEMBL)]
ad_related_con_rosmapP_inc=find_all_connections_multiple_genes(rosmap_trans_rosmap_proteom_inc, rosmap_p_adrelated_inc,2)

write_all_dataframes_to_excel(ad_related_con_rosmapP_inc, 'rosmap_p_inc_v2')



rosmap_t_adrelated_inc=rosmap_trans_inc_unique[which(rosmap_trans_inc_unique %in% AD_related_gene_list$ENSEMBL)]
ad_related_con_rosmapT_inc=find_all_connections_multiple_genes(rosmap_trans_rosmap_proteom_inc, rosmap_t_adrelated_inc,2)

write_all_dataframes_to_excel(ad_related_con_rosmapT_inc, 'rosmap_p_inc_v2')


#rosmap_dec ile APOE etkileşimleri

direct_APOE_connection =find_all_connections(rosmap_trans_dec, 'ENSG00000130203', 1)
writexl::write_xlsx(direct_APOE_connection, 'direct_APOE_connection.xlsx')

second_con_APOE = find_all_connections(rosmap_trans_dec, 'ENSG00000130203', 2)
writexl::write_xlsx(second_con_APOE, 'second_con_APOE.xlsx')

rosmap_trans_all_pattern = rbind(rosmap_trans_dec, rosmap_trans_inc)

second_con_APOE_all_pattern = find_all_connections(rosmap_trans_all_pattern, 'ENSG00000130203', 2)

second_con_APOE_all_pattern_ic_network = find_all_connections_multiple_genes(rosmap_trans_all_pattern, c('ENSG00000130203'), 2)
write_all_dataframes_to_excel(second_con_APOE_all_pattern_ic_network, 'second_con_APOE_ic_network_dahil')


#ROSMAP Transta  APOE ile 2. derece etkileşim ve birden fazla 1.derecedenle etkileşimde olan (sarı olan nodelar)

sari_node_list = c('ENSG00000100399',
                   'ENSG00000101439',
                   'ENSG00000109610',
                   'ENSG00000167614',
                   'ENSG00000171262',
                   'ENSG00000108443',
                   'ENSG00000136986',
                   'ENSG00000107816',
                   'ENSG00000029363',
                   'ENSG00000107560',
                   'ENSG00000180900')

sari_node_APOE_all_pattern = find_all_connections_multiple_genes(rosmap_trans_all_pattern, sari_node_list, 1)
write_all_dataframes_to_excel(sari_node_APOE_all_pattern, 'sari_node_firstDegree_APOE')

sari_node_APOE_all_pattern = find_all_connections_multiple_genes(rosmap_trans_all_pattern, sari_node_list, 2)
write_all_dataframes_to_excel(sari_node_APOE_all_pattern, 'sari_node_secondDegree_APOE')



# coding gene analizi? Kaç tanesis codin gene
coding_gene_list = read.delim("gene_with_protein_product.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
coding_gene_list = coding_gene_list$ensembl_gene_id

long_non_coding = read.delim('RNA_long_non-coding.txt', header = T, sep = '\t',stringsAsFactors = F)
long_non_coding = long_non_coding$ensembl_gene_id


#coding gene list
sum(rosmap_trans_dec_network_unique %in% coding_gene_list)
sum(rosmap_trans_inc_unique %in% coding_gene_list)

sum(rosmap_proteom_dec_unique %in% coding_gene_list)
sum(rosmap_proteom_inc_unique %in% coding_gene_list)

sum(msbb_trans_dec_network_unique %in% coding_gene_list)

second_con_APOE_unique = unique(second_con_APOE$gene1,second_con_APOE$gene2)
sum(second_con_APOE_unique %in% coding_gene_list)
setdiff(second_con_APOE_unique, coding_gene_list)

sari_node_first_degree_unique = c(
        "ENSG00000276075", "ENSG00000268205", "ENSG00000267890", "ENSG00000260916",
        "ENSG00000233927", "ENSG00000232119", "ENSG00000221184", "ENSG00000214160",
        "ENSG00000213015", "ENSG00000205155", "ENSG00000203950", "ENSG00000198380",
        "ENSG00000197860", "ENSG00000189308", "ENSG00000186312", "ENSG00000180900",
        "ENSG00000179292", "ENSG00000179119", "ENSG00000178057", "ENSG00000177683",
        "ENSG00000176340", "ENSG00000171262", "ENSG00000171161", "ENSG00000170906",
        "ENSG00000168569", "ENSG00000167614", "ENSG00000166595", "ENSG00000166471",
        "ENSG00000161048", "ENSG00000161016", "ENSG00000160131", "ENSG00000157500",
        "ENSG00000146232", "ENSG00000146066", "ENSG00000139597", "ENSG00000136986",
        "ENSG00000134152", "ENSG00000133393", "ENSG00000133114", "ENSG00000130255",
        "ENSG00000125534", "ENSG00000122417", "ENSG00000116918", "ENSG00000111911",
        "ENSG00000111237", "ENSG00000109610", "ENSG00000109133", "ENSG00000108443",
        "ENSG00000108107", "ENSG00000108010", "ENSG00000107816", "ENSG00000107560",
        "ENSG00000107317", "ENSG00000107223", "ENSG00000105409", "ENSG00000105254",
        "ENSG00000104154", "ENSG00000103175", "ENSG00000101439", "ENSG00000100399",
        "ENSG00000087258", "ENSG00000077254", "ENSG00000061337", "ENSG00000029363",
        "ENSG00000021574"
)

setdiff(second_con_APOE_unique,coding_gene_list)


#long non coding gene list 
sum(rosmap_trans_dec_network_unique %in% long_non_coding)
sum(rosmap_trans_inc_unique %in% long_non_coding)


rosmap_trans_dec_network_long_non_coding = intersect(rosmap_trans_dec_network_unique, long_non_coding)
rosmap_trans_inc_network_long_non_coding = intersect(rosmap_trans_inc_unique, long_non_coding)

rosmap_trans_all_network_long_non_coding = unique(c(rosmap_trans_dec_network_long_non_coding,rosmap_trans_inc_network_long_non_coding))


rosmap_long_non_coding_1st_degree_con=find_all_connections_combined(rosmap_trans_all_pattern,rosmap_trans_all_network_long_non_coding,1)  
write_all_dataframes_to_excel(rosmap_long_non_coding_1st_degree_con, 'rosmap_long_non_cod_1st_degree_con')


longnon_gene='ENSG00000271119'
rosmap_long_non_coding_1st_degree_con$edges[rosmap_long_non_coding_1st_degree_con$edges$gene1 == longnon_gene|rosmap_long_non_coding_1st_degree_con$edges$gene2 == longnon_gene,]

dd=find_all_connections_multiple_genes(rosmap_trans_all_pattern,rosmap_trans_all_network_long_non_coding,1)

intersect(AD_related_gene_list$ENSEMBL,rosmap_long_non_coding_1st_degree_con$gene_info$gene)


sum(rosmap_proteom_dec_unique %in% long_non_coding)
sum(rosmap_proteom_inc_unique %in% long_non_coding)

sum(msbb_trans_dec_network_unique %in% long_non_coding)

msbb_long_noncoding = intersect(msbb_trans_dec_network_unique, long_non_coding)

msbb_rosmap_longnoncommon=intersect(msbb_long_noncoding, rosmap_trans_dec_network_long_non_coding)






sari_node_first_degree_unique[which(sari_node_first_degree_unique %in% long_non_coding)] #ENSG00000268205


# protein coding  fisher exact test  örnek

# Azalan desen için 2x2 tablo
decreasing_table <- matrix(c(2972, 421,
                             12979, 9502),
                           nrow = 2, byrow = TRUE)
colnames(decreasing_table) <- c("ProteinCoding", "NonCoding")
rownames(decreasing_table) <- c("Decreasing", "NotDecreasing")


lnc_table <- matrix(c(108, 3285,
                      2571, 19910),
                    nrow = 2, byrow = TRUE)
colnames(lnc_table) <- c("LongNonCoding", "OtherGenes")
rownames(lnc_table) <- c("InPattern", "NotInPattern")

lnc_table2 <- matrix(c(36, 466,
                       2643, 22729),
                     nrow = 2, byrow = TRUE)
colnames(lnc_table2) <- c("LongNonCoding", "OtherGenes")
rownames(lnc_table2) <- c("InPattern", "NotInPattern")


fisher.test(lnc_table)
fisher.test(lnc_table)
fisher.test(decreasing_table)


# Tabloyu görüntüle
print(decreasing_table)

# Fisher'ın Kesin Testi
fisher_decreasing <- fisher.test(decreasing_table)
print(fisher_decreasing)







#grafik


rosmap_dec_cor=abs(rosmap_trans_dec$ref_corr)
rosmap_inc_cor=abs(rosmap_trans_inc$ref_corr)

msbb_trans_dec_cor= abs(msbb_trans_dec$ref_corr)
#msbb_trans_inc_cor= abs(msbb_trans_inc$ref_corr)

rosmap_proteom_dec_cor= abs(rosmap_proteom_dec$ref_corr)
rosmap_proteom_inc_cor= abs(rosmap_proteom_inc$ref_corr)


percent_above_99 <- sum(rosmap_dec_cor > 0.99) / length(rosmap_dec_cor) * 100

ggplot(data.frame(abs_corr = rosmap_proteom_dec_cor), aes(x = abs_corr)) +
        geom_histogram(aes(y = (..count..)/sum(..count..) * 100), binwidth = 0.005, fill = "green", color = "black", alpha = 0.7) +
        scale_y_continuous(labels = scales::percent_format(scale = 1)) + # Yüzde işareti eklemek için
        labs(title = "Mutlak Korelasyonların Yüzdesel Dağılımı-Rosmap-p-dec", 
             x = "Mutlak Korelasyon", 
             y = "Yüzde (%)") +
        #xlim(0, max(abs_corr)) + # X eksenini en büyük değere kadar göster
        theme_minimal()


ggplot(data.frame(abs_corr = rosmap_dec_cor), aes(x = abs_corr)) +
        geom_histogram(aes(y = (..count..)/sum(..count..) * 100), binwidth = 0.01, fill = "skyblue", color = "black", alpha = 0.7) +
        scale_y_continuous(labels = scales::percent_format(scale = 1)) + # Yüzde işareti eklemek için
        labs(title = "Mutlak Korelasyonların Yüzdesel Dağılımı-RosmapDec", 
             x = "Mutlak Korelasyon", 
             y = "Yüzde (%)") +
       # xlim(0, max(abs_corr)) + # X eksenini en büyük değere kadar göster
        geom_hline(yintercept = percent_above_99, color = "red", linetype = "dashed", size = 0.5) + # 0.99 çizgisi
        annotate("text", x = 1, y = 1, label = paste0(round(percent_above_99, 2), "%"), color = "red", hjust = 0.1) +
        theme_minimal()





find_all_connections=function(df, gene, degree_limit = Inf) {
        # Initialize the list of found connections
        all_connections = data.frame(gene1 = character(), gene2 = character(), stringsAsFactors = FALSE)
        # Start with the initial gene
        current_genes = c(gene)
        degree = 1
        seen_genes = c(gene)
        
        while (length(current_genes) > 0 & degree <= degree_limit) {
                # Find all pairs involving the current set of genes
                new_connections = df[df$gene1 %in% current_genes | df$gene2 %in% current_genes, ]
                
                # Add new connections to the list of all connections
                all_connections = unique(rbind(all_connections, new_connections))
                
                # Find all genes that have been connected in this iteration
                next_genes = unique(c(new_connections$gene1, new_connections$gene2))
                
                # Remove genes that have already been considered
                current_genes = setdiff(next_genes, seen_genes)
                
                # Add the newly discovered genes to the seen_genes
                seen_genes = unique(c(seen_genes, current_genes))
                
                # Increment degree
                degree = degree + 1
        }
        
        return(all_connections)
}



find_all_connections_multiple_genes=function(df, gene_list, degree_limit = Inf) {
        # Initialize an empty list to store the results for each gene
        all_results = list()
        
        # Loop over each gene in the provided gene list
        for (gene in gene_list) {
                # Initialize the list of found connections for this gene
                all_connections = data.frame(gene1 = character(), gene2 = character(), stringsAsFactors = FALSE)
                # Start with the initial gene
                current_genes = c(gene)
                degree = 1
                seen_genes = c(gene)
                
                while (length(current_genes) > 0 & degree <= degree_limit) {
                        # Find all pairs involving the current set of genes
                        new_connections = df[df$gene1 %in% current_genes | df$gene2 %in% current_genes, ]
                        
                        # Add new connections to the list of all connections
                        all_connections = unique(rbind(all_connections, new_connections))
                        
                        # Find all genes that have been connected in this iteration
                        next_genes = unique(c(new_connections$gene1, new_connections$gene2))
                        
                        # Remove genes that have already been considered
                        current_genes = setdiff(next_genes, seen_genes)
                        
                        # Add the newly discovered genes to the seen_genes
                        seen_genes = unique(c(seen_genes, current_genes))
                        
                        # Increment degree
                        degree = degree + 1
                }
                
                # Store the result for the current gene in the results list
                all_results[[gene]] = all_connections
        }
        
        return(all_results)
}


find_all_connections_multiple_genes = function(df, gene_list, degree_limit = Inf) {
        # Initialize an empty list to store results for each gene
        all_results = list()
        
        # Loop over each gene in the gene list
        for (gene in gene_list) {
                
                # For BFS/bookkeeping
                current_genes = c(gene)       # the frontier for the current "wave"
                seen_genes    = c(gene)       # all genes discovered so far
                degree        = 1
                
                # Expand outward from 'gene', up to 'degree_limit' times
                while (length(current_genes) > 0 & degree <= degree_limit) {
                        
                        # Find all edges where either gene1 or gene2 is in 'current_genes'
                        new_connections = df[
                                df$gene1 %in% current_genes | df$gene2 %in% current_genes, 
                        ]
                        
                        # Identify newly discovered genes from these edges
                        next_genes = unique(c(new_connections$gene1, new_connections$gene2))
                        
                        # 'current_genes' for the next iteration = those we haven't seen yet
                        current_genes = setdiff(next_genes, seen_genes)
                        
                        # Add newly discovered genes to the "seen" list
                        seen_genes = unique(c(seen_genes, current_genes))
                        
                        # Move to the next degree
                        degree = degree + 1
                }
                
                # -- At this point, 'seen_genes' contains the seed gene plus all neighbors 
                #    up through the requested degree. But we still need to pick up edges
                #    *among* these genes themselves. --
                
                # Filter the original df to keep only edges between genes in 'seen_genes'
                final_edges = df[
                        df$gene1 %in% seen_genes & df$gene2 %in% seen_genes, 
                ]
                
                # Store in a list keyed by the seed gene
                all_results[[gene]] = final_edges
        }
        
        return(all_results)
}

find_all_connections_combined <- function(df, gene_list, degree_limit = Inf) {
        # Kenarları ve genleri toplamak için başlangıç değişkenleri
        all_edges <- data.frame(gene1 = character(), gene2 = character(), stringsAsFactors = FALSE)
        all_nodes <- character()
        
        # Her gen için BFS ile genişleme yapıyoruz
        for (gene in gene_list) {
                current_genes <- c(gene)  # o anki dalga (wave) genleri
                seen_genes <- c(gene)     # daha önce görülmüş genler
                degree <- 1
                
                while (length(current_genes) > 0 && degree <= degree_limit) {
                        # current_genes içeren kenarları seç
                        new_connections <- df[df$gene1 %in% current_genes | df$gene2 %in% current_genes, ]
                        # Global kenar listesine ekle (duplicate’ları engellemek için unique)
                        all_edges <- unique(rbind(all_edges, new_connections))
                        
                        # Bu turda bulunan tüm genleri belirle
                        next_genes <- unique(c(new_connections$gene1, new_connections$gene2))
                        # Daha önce görülmeyenleri al: sonraki dalga
                        current_genes <- setdiff(next_genes, seen_genes)
                        # Görülen genleri güncelle
                        seen_genes <- unique(c(seen_genes, current_genes))
                        
                        # Dereceyi arttır
                        degree <- degree + 1
                }
                
                # Bulunan genleri global listeye ekle
                all_nodes <- union(all_nodes, seen_genes)
        }
        
        # Tüm bulunan genler arasındaki kenarları orijinal df üzerinden alalım 
        combined_edges <- df[
                df$gene1 %in% all_nodes & df$gene2 %in% all_nodes, 
        ]
        
        # Ana genlerin (gene_list) mutlaka dahil olduğundan emin olalım
        all_nodes <- union(all_nodes, gene_list)
        
        # gene_info data.frame’i: 1. kolon tüm bulunan gen, 2. kolon main/side bilgisi
        gene_info <- data.frame(gene = all_nodes, stringsAsFactors = FALSE)
        gene_info$type <- ifelse(gene_info$gene %in% gene_list, "main", "side")
        
        # İki tabloyu bir liste olarak döndürüyoruz
        return(list(edges = combined_edges, gene_info = gene_info))
}




write_all_dataframes_to_excel=function(data_list, base_path, file_extension = ".xlsx") {
        # Check if the directory exists; if not, create it
        if (!dir.exists(base_path)) {
                dir.create(base_path, recursive = TRUE)
                cat("Directory created:", base_path, "\n")
        }
        
        # Iterate over each element in the list
        for (element_name in names(data_list)) {
                # Create the file name based on element_name
                file_name <- paste0(element_name, file_extension)
                
                # Construct the full file path
                file_path <- file.path(base_path, file_name)
                
                # Save the DataFrame as an xlsx file
                writexl::write_xlsx(data_list[[element_name]], file_path)
                
                # Print the saved file path for verification
                cat("Saved:", file_path, "\n")
        }
}

