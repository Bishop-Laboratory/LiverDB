# This R script requires:
#
# - GTF_FILE refers to a valid .gtf file
# - COUNTS_DIR directory populated with raw counts files for all samples
#   listed in METADATA_CSV
#
# A successful run of the Nextflow script "pipeline.nf" will generate these
# requirements.

#-------------------------------------------------------------------------------
# LIBRARIES
#-------------------------------------------------------------------------------
library(edgeR)
library(tidyverse)

#-------------------------------------------------------------------------------
# CONSTANTS
#-------------------------------------------------------------------------------

COUNTS_DIR <- "raw_counts"
COUNTS_SKIP_LINES <- 4

GTF_FILE <- "refs/Homo_sapiens.GRCh38.103.gtf"

METADATA_CSV <- "metadata/metadata.csv"
CONTRASTS_CSV <- "metadata/contrasts.csv"

EXP_CSV_SUFFIX <- "_gene_exp.csv"
DEG_CSV_SUFFIX <- "_degs.csv"


#-------------------------------------------------------------------------------
# HELPER FUNCTIONS
#-------------------------------------------------------------------------------

mat_longify <- function(x) {
  as.data.frame(x) %>%
    rownames_to_column("gene_id") %>%
    pivot_longer(cols = contains("SRR"))
}

# http://luisvalesilva.com/datasimple/rna-seq_units.html
rpkm_to_tpm <- function(x){
  rpkm.sum <- colSums(x)
  return(t(t(x) / (1e-06 * rpkm.sum)))
}


#-------------------------------------------------------------------------------
# MAIN SCRIPT
#-------------------------------------------------------------------------------

all_counts <- list.files(COUNTS_DIR)
metadata <- read_csv(METADATA_CSV)
all_contrasts <- read_csv(CONTRASTS_CSV)

# Get gene lengths
txdb <- GenomicFeatures::makeTxDbFromGFF(file = GTF_FILE)
gene_lengths <- GenomicFeatures::transcriptsBy(txdb, "gene") %>%
  GenomicRanges::reduce() %>%
  GenomicRanges::width() %>%
  sum()

# For each unique study, generate expression level and DEG .csv files
unique_studies <- metadata %>% distinct(study_id) %>% pull()
for (study_id in unique_studies) {
  # Filter samples relevant to current study
  samples <- metadata %>% filter(study_id == !!study_id) %>% pull(sample_id)
  study_counts <- all_counts[str_extract(all_counts, "SRR\\d+") %in% samples]
  
  # Build raw count matrix
  mat <- lapply(study_counts, function(x) {
    sample_id <- str_extract(x, "SRR\\d+")
    strand <- metadata$stranded[metadata$sample_id == sample_id]
    count_path <- file.path(COUNTS_DIR, x)
    
    read_tsv(count_path,
      skip = COUNTS_SKIP_LINES, 
      col_names = c("gene_id", "unstranded", "forward", "reverse")
    ) %>%
      select(gene_id, contains(!!strand)) %>%
      rename(counts = contains(!!strand)) %>%
      mutate(sample_id = !!sample_id)
  }) %>%
    bind_rows() %>%
    pivot_wider(names_from = sample_id, values_from = counts) %>%
    column_to_rownames(var="gene_id") %>%
    as.matrix()
  
  # Compute CPM, RPKM and TPM
  cpms <- mat %>%
    DGEList() %>%
    calcNormFactors() %>%
    cpm() %>% 
    mat_longify() %>%
    rename(sample_id = name, cpm = value)

  rpkms_mat <- rpkm(mat, gene_lengths[rownames(mat)])
  rpkms <- rpkms_mat %>%
    mat_longify() %>%
    rename(sample_id = name, rpkm = value)

  tpms <- tpm_from_rpkm(rpkms_mat) %>%
    mat_longify() %>%
    rename(sample_id = name, tpm = value)
  
  # Write CPM, RPKM and TPM to .csv
  cpms %>% 
    inner_join(rpkms, by = c("gene_id", "sample_id")) %>% 
    inner_join(tpms, by = c("gene_id", "sample_id")) %>% 
    write_csv(paste0(study_id, EXP_CSV_SUFFIX))
  
  # Remove large unneeded R objects and free up some memory
  rm(cpms, rpkms, rpkms_mat, tpms)
  gc()
  
  
  # Get contrasts relevant to current study
  study_contrasts <- all_contrasts %>%
    filter(study_id == !!study_id)
  
  # Get DEGs for each contrast and write to .csv
  apply(study_contrasts, 1, function(x) {
    # Define current contrast
    numerator <-  x[["numerator"]]
    denominator <- x[["denominator"]]
    
    # Get relevant samples according to current contrast
    deg_meta <- metadata %>%
      filter(study_id == !!study_id) %>%
      mutate(group = case_when(
        str_starts(condition, numerator) ~ 1,
        str_starts(condition, denominator) ~ 2,
        TRUE ~ 0
      )) %>%
      filter(group != 0)
    
    groups <- factor(deg_meta$group)
    design <- model.matrix(~0 + groups)
    
    # DGE analysis
    mat[, deg_meta$sample_id] %>%
      DGEList(group = groups) %>%  
      calcNormFactors() %>%
      estimateGLMCommonDisp(design) %>%
      estimateGLMTrendedDisp(design) %>%
      estimateGLMTagwiseDisp(design) %>%
      glmFit(design) %>%
      glmLRT(contrast = c(-1, 1)) %>%
      topTags(n="Inf") %>% 
      pluck("table") %>%
      mutate(
        numerator = numerator,
        denominator, denominator
      ) %>% 
      rownames_to_column("gene_id") %>% 
      select(gene_id, numerator, denominator, logFC, FDR)
  }) %>% 
    bind_rows() %>% 
    write_csv(paste0(study_id, DEG_CSV_SUFFIX))
  
  # Remove large unneeded R objects and free up some memory
  rm(mat)
  gc()

}

