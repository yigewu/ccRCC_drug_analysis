# Yige Wu @ WashU 2021 Mar

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input average protein
avg_pro_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/avg_protein_by_model/20210303.v1/Average_Protein.By_Model.20210303.v1.tsv")

# expand by gene and filter -----------------------------------------------
## filter by NA
avg_pro_filtered_df <- avg_pro_df %>%
  filter(!is.na(RESL5) & !is.na(RESL10))
## preprocess
avg_pro_filtered_df[avg_pro_filtered_df$PG.Genes == "H1-0", "PG.Genes"] <- "H1-0;H1-0"
## get repetition index
idx_rep <- sapply(avg_pro_filtered_df$PG.Genes, function(gene_string) {
  genes_vec <- str_split(string = gene_string, pattern = ";")[[1]]
  len_rep <- length(genes_vec)
  return(len_rep)
})
genes_uniq_vec <- sapply(avg_pro_filtered_df$PG.Genes, function(gene_string) {
  genes_vec <- str_split(string = gene_string, pattern = ";")[[1]]
  return(genes_vec)
})
proteinnames_uniq_vc <- sapply(avg_pro_filtered_df$PG.ProteinNames, function(gene_string) {
  genes_vec <- str_split(string = gene_string, pattern = ";")[[1]]
  return(genes_vec)
})
proteinaccessions_uniq_vc <- sapply(avg_pro_filtered_df$PG.ProteinAccessions, function(gene_string) {
  genes_vec <- str_split(string = gene_string, pattern = ";")[[1]]
  return(genes_vec)
})

## create new data frame
avg_pro_new_df <- avg_pro_filtered_df[rep(1:nrow(avg_pro_filtered_df), idx_rep),]
avg_pro_new_df$PG.Gene <- unlist(genes_uniq_vec)
avg_pro_new_df$PG.ProteinName <- unlist(proteinnames_uniq_vc)
avg_pro_new_df$PG.ProteinAccession <- unlist(proteinaccessions_uniq_vc)
## select columns
avg_pro_new_df <- avg_pro_new_df %>%
  select(PG.Gene, PG.ProteinName, PG.ProteinAccession, PG.ProteinNames, RESL10, RESL5) %>%
  mutate(log2Intensity_FC_RESL10vs5 = (RESL10 - RESL5)) %>%
  mutate(RESL10_vs_RESL5.protein = ifelse(log2Intensity_FC_RESL10vs5 > 0, "Up", "Down"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "RESL10_vs_RESL5.", "Protein.", run_id, ".tsv")
write.table(x = avg_pro_new_df, file = file2write, quote = F, sep = "\t", row.names = F)
