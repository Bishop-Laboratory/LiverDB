# This R script requires:
#
# - getwd() is the same as the location of this script
# - GTF_FILE refers to a valid .gtf file
# - COUNTS_DIR directory populated with raw counts files for all samples
#   listed in METADATA_CSV
# 
# A successful run of the Nextflow script "pipeline.nf" will generate the file
# requirements above.

#-------------------------------------------------------------------------------
# LIBRARIES
#-------------------------------------------------------------------------------
library(edgeR)
library(tidyverse)

#-------------------------------------------------------------------------------
# CONSTANTS
#-------------------------------------------------------------------------------

COUNTS_DIR <- "raw_counts"
COUNTS_SUFFIX <- "_ReadsPerGene.out.tab"
COUNTS_SKIP_LINES <- 4

GTF_FILE <- "refs/Homo_sapiens.GRCh38.103.gtf"

METADATA_CSV <- "metadata/metadata.csv"
CONTRASTS_CSV <- "metadata/contrasts.csv"

EXP_OUTFILE_SUFFIX <- "_gene_exp.csv.gz"
DEG_OUTFILE_SUFFIX <- "_degs.csv.gz"


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

metadata <- read_csv(METADATA_CSV)
all_contrasts <- read_csv(CONTRASTS_CSV)

# Get gene lengths
txdb <- GenomicFeatures::makeTxDbFromGFF(file = GTF_FILE)
gene_lengths <- GenomicFeatures::transcriptsBy(txdb, "gene") %>%
  GenomicRanges::reduce() %>%
  GenomicRanges::width() %>%
  sum()

# For each unique study, generate .csv.gz files for CPM/RKPM/TPM and DEGs
unique_studies <- metadata %>% distinct(study_id) %>% pull()
for (study_id in unique_studies) {
  exp_filename <- paste0(study_id, EXP_OUTFILE_SUFFIX)
  deg_filename <- paste0(study_id, DEG_OUTFILE_SUFFIX)
  
  # Filter sample IDs list relevant to current study
  samples <- metadata %>% filter(study_id == !!study_id) %>% pull(sample_id)
  
  # Build raw count matrix
  mat <- lapply(samples, function(id) {
    strand <- metadata %>% filter(sample_id == id) %>% pull(stranded)
    count_path <- file.path(COUNTS_DIR, paste0(id, COUNTS_SUFFIX))
    
    read_tsv(count_path,
      skip = COUNTS_SKIP_LINES, 
      col_names = c("gene_id", "unstranded", "forward", "reverse")
    ) %>%
      select(gene_id, contains(!!strand)) %>%
      rename(counts = contains(!!strand)) %>%
      mutate(sample_id = id)
  }) %>%
    bind_rows() %>%
    pivot_wider(names_from = sample_id, values_from = counts) %>%
    column_to_rownames("gene_id") %>%
    as.matrix()
  
  # Compute CPM, RPKM, TPM
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

  tpms <- rpkm_to_tpm(rpkms_mat) %>%
    mat_longify() %>%
    rename(sample_id = name, tpm = value)
  
  # Write CPM, RPKM and TPM to .csv.gz
  cpms %>% 
    inner_join(rpkms, by = c("gene_id", "sample_id")) %>% 
    inner_join(tpms, by = c("gene_id", "sample_id")) %>% 
    write_csv(exp_filename)
  
  # Remove large unneeded R objects and free up some memory
  rm(cpms, rpkms, rpkms_mat, tpms)
  gc()
  
  
  # Get contrasts relevant to current study
  study_contrasts <- all_contrasts %>%
    filter(study_id == !!study_id)
  
  # Get DEGs for each contrast and write to .csv.gz
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
        denominator = denominator
      ) %>% 
      rownames_to_column("gene_id") %>% 
      select(gene_id, numerator, denominator, logFC, FDR)
  }) %>% 
    bind_rows() %>% 
    write_csv(deg_filename)
  
  # Remove large unneeded R objects and free up some memory
  rm(mat)
  gc()

}

