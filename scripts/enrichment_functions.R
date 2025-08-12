# File: enrichment_functions.R
# This script contains functions for performing Gene Ontology (GO) and KEGG pathway
# enrichment analysis using the clusterProfiler package.

# 1. Load Libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)

# 2. Enrichment Analysis Functions

# Function to perform Gene Ontology (GO) enrichment analysis.
# This function is based on the user's `enrichment` function.
enrichment_go <- function(gene_list, ontology_type, title) {
  # Perform GO enrichment analysis using clusterProfiler::enrichGO
  enrichment_result <- enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = 'ENSEMBL',
    ont = ontology_type,
    pvalueCutoff = 0.05,
    pAdjustMethod = 'BH',
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE,
    pool = FALSE
  )
  return(enrichment_result)
}

# Function to perform KEGG pathway enrichment analysis.
# This function is based on the user's `enrichment_KEGG` function.
enrichment_kegg <- function(gene_list, title) {
  # Convert ENSEMBL IDs to Entrez IDs, which is required for KEGG analysis
  entrezID <- unlist(mget(gene_list, envir = org.Hs.egENSEMBL2EG, ifnotfound = NA))
  entrezID <- na.omit(entrezID)
  
  # Perform KEGG enrichment analysis
  kegg_result <- enrichKEGG(
    gene = entrezID,
    pvalueCutoff = 0.05
  )
  return(kegg_result)
}

# Function to visualize enrichment results.
# This function is derived from the plotting logic in `patternYeni.R`.
plot_enrichment_results <- function(combined_data) {
  # Add -log10(pvalue) for sizing the points.
  combined_data$neg_log10_pvalue <- -log10(combined_data$pvalue)
  
  # Set the order of patterns for the plot.
  combined_data$Pattern <- factor(
    combined_data$Pattern,
    levels = c("Decreasing Pattern", "Increasing Pattern")
  )
  
  # Define a specific order for the descriptions to maintain consistency across plots.
  my_order <- c(
    "Alzheimer disease",
    "aerobic respiration",
    "sphingomyelin biosynthetic process",
    "vesicle organization",
    "neuron to neuron synapse", "asymmetric synapse", 'postsynaptic density', 'Endocytosis',
    'Ubiquitin mediated proteolysis', 'protein polyubiquitination',
    'dendrite development', 'transport along microtubule', 'macroautophagy',
    'GTPase activity', 'Oxidative phosphorylation', 'respiratory chain complex',
    'purine nucleoside triphosphate biosynthetic process', 'ATP biosynthetic process',
    'NADH dehydrogenase (ubiquinone) activity', 'proton motive force-driven mitochondrial ATP synthesis',
    'mitochondrial protein-containing complex', 'mitochondrial ATP synthesis coupled electron transport'
  )
  
  # Reorder descriptions based on the defined order.
  combined_data$Description <- factor(
    combined_data$Description,
    levels = rev(my_order)
  )
  
  # Create the dot plot.
  p <- ggplot(combined_data, aes(x = Pattern, y = Description)) +
    geom_point(aes(size = neg_log10_pvalue), colour = "firebrick") +
    facet_grid(. ~ Pattern, scales = "free_x", space = "free_x") +
    scale_size_continuous(name = "-log10(p-value)", range = c(3, 8)) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "lightgrey", colour = "black"),
      strip.text = element_text(face = "bold", size = 12),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(1, "lines"),
      legend.position = "right"
    )
  
  return(p)
}