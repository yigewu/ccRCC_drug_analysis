# Yige Wu @ WashU 2020 May

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
## input mutations for RCC PDX
# mut_df <- fread(data.table = F, input = "./Resources/Analysis_Results/mutation/generate_bulk_mutation_table/20200526.v1/RCC_PDX.somaticMut.20200526.v1.tsv")
mut_df <- fread(data.table = F, input = "./Resources/Analysis_Results/mutation/generate_bulk_mutation_table/20200528.v1/RCC_PDX.somaticMut.20200528.v1.tsv")

# extract SMGs and transform into wide data frame----------------------------------------------
mut_long_df <- mut_df %>%
  filter(Hugo_Symbol %in% ccRCC_SMGs)
mut_long_unique_df <- mut_long_df %>%
  select(Analysis_ID, Hugo_Symbol, Variant_Classification) %>%
  unique()
mut_wide_df <- dcast(data = mut_long_unique_df, formula = Analysis_ID ~ Hugo_Symbol, value.var = "Variant_Classification")
mut_wide_df[is.na(mut_wide_df)] <- "None"
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "RCC_PDX.", "somaticMut.", "SMG.", "Wide.", run_id, ".tsv")
write.table(x = mut_wide_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "RCC_PDX.", "somaticMut.", "SMG.", "Long.", run_id, ".tsv")
write.table(x = mut_long_df, file = file2write, quote = F, sep = "\t", row.names = F)
