# Yige Wu @ WashU 2021 Jan

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
exp_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/preprocess/normalize_DIA_protein_data/20210111.v1/RCC_PDX.DIA_Protein.Log2.QuantileNormalized.20210111.v1.tsv")
## input the sample info
meta_data_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")

# expand gene symbols -----------------------------------------------
## preprocess
exp_df[exp_df$PG.Genes == "H1-0", "PG.Genes"] <- "H1-0;H1-0"
## get repetition index
idx_rep <- sapply(exp_df$PG.Genes, function(gene_string) {
  genes_vec <- str_split(string = gene_string, pattern = ";")[[1]]
  len_rep <- length(genes_vec)
  return(len_rep)
})
genes_uniq_vec <- sapply(exp_df$PG.Genes, function(gene_string) {
  genes_vec <- str_split(string = gene_string, pattern = ";")[[1]]
  return(genes_vec)
})
proteinnames_uniq_vc <- sapply(exp_df$PG.ProteinNames, function(gene_string) {
  genes_vec <- str_split(string = gene_string, pattern = ";")[[1]]
  return(genes_vec)
})
proteinaccessions_uniq_vc <- sapply(exp_df$PG.ProteinAccessions, function(gene_string) {
  genes_vec <- str_split(string = gene_string, pattern = ";")[[1]]
  return(genes_vec)
})

## create new data frame
exp_process_df <- exp_df[rep(1:nrow(exp_df), idx_rep),]
exp_process_df$PG.Gene <- unlist(genes_uniq_vec)
exp_process_df$PG.ProteinName <- unlist(proteinnames_uniq_vc)
exp_process_df$PG.ProteinAccession <- unlist(proteinaccessions_uniq_vc)

# process -------------------------------------------------------------
exp_idx_df <- exp_process_df[, c("PG.Gene", "PG.ProteinName", "PG.ProteinDescriptions", "PG.ProteinAccession")]
exp_data_df <- exp_process_df[, meta_data_df$`Sample ID`]
exp_zscore_t_mat <- apply(exp_data_df, 1, FUN = function(x) {
  (x-median(x, na.rm = T)) / sd(x, na.rm = T)
})
exp_zscore_df <- data.frame(t(exp_zscore_t_mat))
exp_zscore_df <- cbind(exp_idx_df, exp_zscore_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Protein_zscore.", run_id, ".tsv")
write.table(x = exp_zscore_df, file = file2write, quote = F, sep = "\t", row.names = F)
