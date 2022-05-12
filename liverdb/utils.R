# UI function to make ? button
helpButton <- function(message) {
  return(
    add_prompt(
      ui_element = span(HTML('<i class="fa fa-question-circle"></i>')),
      message = message, position = "right"
    )
  )
}

#' Make headers
makeHeaders <- function(title, message, fs=1.3) {
  tagList(
    span(span(title, style=paste0("font-size: ", fs, "em;")), helpButton(message))
  )
}


#' Makes the global data for the app
# makeGlobalData <- function(APP_DATA) {
#   exp <- readr::read_csv("https://fibrodb-data.s3.amazonaws.com/gene_exp.csv.gz")
#   samples <- readr::read_csv("https://fibrodb-data.s3.amazonaws.com/samples.csv")
#   contrasts <- readr::read_csv("https://fibrodb-data.s3.amazonaws.com/contrasts.csv")
#   degs <- readr::read_csv("https://fibrodb-data.s3.amazonaws.com/degs.csv.gz")
#   app_data <- list(
#     exp=exp, samples=samples, contrasts=contrasts, degs=degs
#   )
#   ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
#   G_list <- biomaRt::getBM(attributes= c("ensembl_gene_id","hgnc_symbol"), mart= ensembl)
#   ens2sym <- G_list %>%
#     filter(hgnc_symbol != "") %>%
#     rename(gene_id=ensembl_gene_id, gene_name=hgnc_symbol)
#   results <- full_join(
#     exp, samples
#   ) %>%
#     full_join(
#       degs
#     ) %>%
#     full_join(
#       contrasts
#     ) %>%
#     inner_join(
#       ens2sym
#     )
#   results <- results %>%
#     relocate(gene_name) %>%
#     arrange(padj)
#   results_show <- results %>%
#     select(gene_name, study_id, numerator,
#            denominator, fc, padj) %>%
#     distinct(gene_name, study_id, .keep_all = TRUE)
#   # Add in pathway analysis
#   eres <- results_show %>%
#     filter(! is.na(padj) & padj < .05) %>%
#     group_by(study_id) %>%
#     {setNames(group_split(.), group_keys(.)[[1]])} %>%
#     lapply(
#       function(x) {
#         resup <- enrichR::enrichr(
#           genes = x$gene_name[x$fc > 1], databases = "KEGG_2019_Human"
#         ) %>%
#           purrr::pluck("KEGG_2019_Human") %>%
#           mutate(group = "Over-expressed")
#         resdn <- enrichR::enrichr(
#           genes = x$gene_name[x$fc < -1], databases = "KEGG_2019_Human"
#         ) %>%
#           purrr::pluck("KEGG_2019_Human") %>%
#           mutate(group = "Under-expressed")
#         bind_rows(resup, resdn)
#       }
#     )
#   saveRDS(eres, file = "eres.rds", compress = "xz")
#   saveRDS(results, file = "app_data.rds", compress = "xz")
#   return(app_data)
# }



## Create eres results link
makeERES <- function() {
  eres <- readRDS("eres.rds")
  results <- readRDS("app_data.rds")
  results <- results %>% 
    mutate(
      fc = case_when(
        study_id == "GSE149413" ~ -1 * fc, TRUE ~ fc
      ),
      numerator = case_when(
        study_id == "GSE149413" ~ "Thrombus", TRUE ~ numerator
      ),
      denominator = case_when(
        study_id == "GSE149413" ~ "Adventitia", TRUE ~ denominator
      )
    )
  eres$GSE149413 <- eres$GSE149413 %>% 
    mutate(
      group = ifelse(group == "Under-expressed", "Over-expressed", "Under-expressed")
    )
  contrasts_new <- results %>% filter(! is.na(padj)) %>% dplyr::select(study_id, numerator, denominator) %>% distinct(study_id, .keep_all = TRUE)
  towrite <- lapply(names(eres), function(nx){
    x <- eres[[nx]]
    x$study_id <- nx
    # Remove "Old" p value columns and change column order
    dplyr::select(x, -contains("Old")) %>% 
      dplyr::relocate(study_id, group)
  }) %>% bind_rows() 
  left_join(contrasts_new, towrite) %>% write_csv(file = "enrichment_res.csv")
  system("aws s3 cp enrichment_res.csv s3://fibrodb-data/")
  file.remove("enrichment_res.csv")
}

