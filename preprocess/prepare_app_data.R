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
library(tidyverse)

#-------------------------------------------------------------------------------
# CONSTANTS
#-------------------------------------------------------------------------------

APP_DATA_FILENAME <- "app_data.rds"

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
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ens2sym <- getBM(
  attributes = c("ensembl_gene_id","hgnc_symbol"),
  mart = ensembl
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
  read_csv(filename) %>% 
    inner_join(ens2sym, by = c("gene_id")) %>% 
    relocate(gene_name) %>% 
    select(-gene_id)
})

# Generate app data
app_data <-  list(
  metadata = metadata_tbl,
  contrasts = contrasts_tbl,
  degs = degs,
  exps = exps
)

saveRDS(app_data, file = APP_DATA_FILENAME)


