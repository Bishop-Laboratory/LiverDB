# This R script requires:
#
# - getwd() is the same as the location of this script
# - For each study, .csv.gz files for DEGs and expression levels exist in the
#   same location of this script
# 
# A successful run of the R script "process_raw_counts.R" and "enrichr_res.R"
# will generate the file requirements above.

#-------------------------------------------------------------------------------
# LIBRARIES
#-------------------------------------------------------------------------------
library(biomaRt)
library(tidyverse)

#-------------------------------------------------------------------------------
# CONSTANTS
#-------------------------------------------------------------------------------

APP_DATA_FILENAME <- "app_data.rds"

ERES_CSV <- "enrichr_res.csv.gz"
METADATA_CSV <- "metadata/metadata.csv"
CONTRASTS_CSV <- "metadata/contrasts.csv"
DEG_SUFFIX <- "_degs.csv.gz"
EXP_SUFFIX <- "_gene_exp.csv.gz"


#-------------------------------------------------------------------------------
# MAIN SCRIPT
#-------------------------------------------------------------------------------

# Get study IDs from contrasts table
metadata_tbl <- read_csv(METADATA_CSV)
contrasts_tbl <- read_csv(CONTRASTS_CSV)
study_ids <- contrasts_tbl %>% 
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

# Get DEGs and expression levels, and add gene symbols with inner join (removes
# gene IDs without a symbol)
degs <- lapply(study_ids, function(id) {
  filename <- paste0(id, DEG_SUFFIX)
  read_csv(filename) %>% 
    inner_join(ens2sym, by = c("gene_id")) %>% 
    relocate(gene_name) %>% 
    select(-gene_id)
})

exps <- lapply(study_ids, function(id) {
  filename <- paste0(id, EXP_SUFFIX)
  
  cpm <- read_csv(filename) %>% 
    select(gene_id, sample_id, cpm) %>% 
    pivot_wider(names_from = sample_id, values_from = cpm) %>% 
    inner_join(ens2sym, by = c("gene_id")) %>% 
    relocate(gene_name) %>% 
    select(-gene_id)
  
  rpkm <- read_csv(filename) %>% 
    select(gene_id, sample_id, rpkm) %>% 
    pivot_wider(names_from = sample_id, values_from = rpkm) %>% 
    inner_join(ens2sym, by = c("gene_id")) %>% 
    relocate(gene_name) %>% 
    select(-gene_id)
  
  tpm <- read_csv(filename) %>% 
    select(gene_id, sample_id, tpm) %>% 
    pivot_wider(names_from = sample_id, values_from = tpm) %>% 
    inner_join(ens2sym, by = c("gene_id")) %>% 
    relocate(gene_name) %>% 
    select(-gene_id)
  
  list(cpm = cpm, rpkm = rpkm, tpm = tpm)
})

# Get enrichr results and convert to list
eres_tbl <- read_csv(ERES_CSV)

study_contrasts <- unique(eres_tbl$Study_Contrast)
names(study_contrasts) <- study_contrasts

eres <- lapply(as.list(study_contrasts), function(x) {
  eres_tbl %>% 
    filter(Study_Contrast == x) %>% 
    select(-Study_Contrast)
})

# Generate app data
app_data <-  list(
  metadata = metadata_tbl,
  contrasts = contrasts_tbl,
  degs = degs,
  exps = exps,
  eres = eres
)

saveRDS(app_data, file = APP_DATA_FILENAME)


