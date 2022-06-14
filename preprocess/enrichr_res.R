# This R script requires:
#
# - getwd() is the same as the location of this script
# - For each study, .csv.gz files for DEGs and expression levels exist in the
#   same location of this script
# 
# A successful run of the R script "process_raw_counts.R" will generate the file
# requirements above.

#-------------------------------------------------------------------------------
# LIBRARIES
#-------------------------------------------------------------------------------
library(biomaRt)
library(enrichR)
library(tidyverse)

#-------------------------------------------------------------------------------
# CONSTANTS
#-------------------------------------------------------------------------------

OUTPUT_FILENAME <- "enrichr_res.csv.gz"

METADATA_CSV <- "metadata/metadata.csv"
DEG_SUFFIX <- "_degs.csv.gz"


#-------------------------------------------------------------------------------
# MAIN SCRIPT
#-------------------------------------------------------------------------------

# Get study IDs
study_ids <- read_csv(METADATA_CSV) %>% 
  pull(study_id) %>% 
  unique()
names(study_ids) <- study_ids

# Get gene symbols from BioMart
ens2sym <- getBM(
  attributes = c("ensembl_gene_id","hgnc_symbol"),
  mart = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
) %>% 
  filter(hgnc_symbol != "") %>%
  rename(gene_id = ensembl_gene_id, gene_name = hgnc_symbol)

# Get DEG tables with gene symbols, and filter and label sig DEGs
degs <- lapply(study_ids, function(study){
  filename <- paste0(study, DEG_SUFFIX)
  read_csv(filename) %>% 
    inner_join(ens2sym, by = c("gene_id")) %>% 
    mutate(study_id = study)
}) %>% 
  bind_rows() %>% 
  filter(!is.na(FDR) & FDR < .01 & abs(logFC) > 1) %>%
  unite("group", c("study_id", "numerator", "denominator"))

# Get unique groups of enrichr runs
unique_groups <- degs %>%
  pull(group) %>%
  unique()
names(unique_groups) <- unique_groups

# Send genes to enrichr
eres <- lapply(unique_groups, function(x) {
  up_genes <- degs %>% 
    filter(group == x & logFC > 1) %>% 
    pull(gene_name)
  dn_genes <- degs %>% 
    filter(group == x & logFC < -1) %>% 
    pull(gene_name)
  
  resup <- enrichr(up_genes, databases = "KEGG_2019_Human") %>%
    pluck("KEGG_2019_Human") %>%
    mutate(group = "Over-expressed")
  resdn <- enrichr(dn_genes, databases = "KEGG_2019_Human") %>%
    pluck("KEGG_2019_Human") %>%
    mutate(group = "Under-expressed")
  
  bind_rows(resup, resdn)
})

lapply(as.list(names(eres)), function(x) {
  eres[[x]] %>% 
    mutate(Study_Contrast = x) %>% 
    relocate(Study_Contrast)
}) %>% 
  bind_rows() %>% 
  write_csv(OUTPUT_FILENAME)


