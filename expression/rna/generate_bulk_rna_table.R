# Yige Wu @ WashU 2020 Feb

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

# input dependencies --------------------------------------------
## input all ccRCC samples data info
data_status_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_status/write_bulk_data_status/20200612.v1/RCC_PDX_Related_Samples.PDXNet_B1_9.HTAN_B1.Data_Status.20200612.v1.tsv")
## input batch 1-8 gene expression
rna_exp_b1_8_df <- fread("./Resources/Bulk_Processed_Data/Data_Files/batch1_8/GeneExp/geneExp.ccrcc.kallisto.tpm.20200316.tsv", data.table = F)
## input batch 9 gene expression
rna_exp_b9_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Data_Files/batch9/geneExp/b9.pdx.kallisto.tpm.gene_level.tsv")

# rename batch 1-8 --------------------------------------------------
## get the colnames that are actual rcc pdx analysis ids
colnames_keep <- colnames(rna_exp_b1_8_df)[colnames(rna_exp_b1_8_df) %in% c("Name", data_status_df$Analysis_ID.RNA[!is.na(data_status_df$Analysis_ID.RNA)])]
rna_exp_b1_8_new_df <- rna_exp_b1_8_df[, colnames_keep]
rna_exp_b1_8_new_df <- rna_exp_b1_8_new_df %>%
  rename(gene_symbol = Name)
nrow(rna_exp_b1_8_new_df)
rna_exp_b1_8_new_df$gene_symbol

# rename b9 ---------------------------------------------------------------
## get the colnames that are actual rcc pdx analysis ids
colnames(rna_exp_b9_df)
colnames_keep <- colnames(rna_exp_b9_df)[colnames(rna_exp_b9_df) %in% c("V1", data_status_df$Analysis_ID.RNA[!is.na(data_status_df$Analysis_ID.RNA)])]
colnames_keep
rna_exp_b9_new_df <- rna_exp_b9_df[, colnames_keep]
rna_exp_b9_new_df <- rna_exp_b9_new_df %>%
  rename(gene_symbol = V1)
nrow(rna_exp_b9_new_df)

# merge -------------------------------------------------------------------
rna_exp_b1_9_df <- merge(rna_exp_b9_new_df, rna_exp_b1_8_new_df, by = c("gene_symbol"), all = T)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "RCC_PDX.", "geneExp.", run_id, ".tsv")
write.table(x = rna_exp_b1_9_df, file = file2write, quote = F, sep = "\t", row.names = F)



