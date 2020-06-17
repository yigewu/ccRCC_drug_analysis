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
## input the CNV data for RCC PDX
genecnv_df <- fread(data.table = F, input = "./Resources/Analysis_Results/copy_number/generate_gene_level_copy_number_table/20200526.v1/RCC_PDX.Gene_Level_CNV.20200526.v1.tsv")
## input known CNV genes
knowncnvgenes_df <- readxl::read_xlsx(path = "./Resources/Known_Genetic_Alterations/Known_CNV.20200505.v1.xlsx", sheet = "Genes")

# extract specific CNV genes and transform into wide data frame----------------------------------------------
genecnv_long_df <- genecnv_df %>%
  filter(gene %in% knowncnvgenes_df$Gene_Symbol)
genecnv_wide_log2_df <- dcast(data = genecnv_long_df, formula = sample ~ gene, value.var = "log2")
genecnv_wide_cn_df <- dcast(data = genecnv_long_df, formula = sample ~ gene, value.var = "cn")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "RCC_PDX.Gene_Level_CNV.", "Known_CNV_Genes.", "Wide.", "Log2.", run_id, ".tsv")
write.table(x = genecnv_wide_log2_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "RCC_PDX.Gene_Level_CNV.", "Known_CNV_Genes.", "Wide.", "CN.", run_id, ".tsv")
write.table(x = genecnv_wide_cn_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "RCC_PDX.Gene_Level_CNV.", "Known_CNV_Genes.", "Long.", run_id, ".tsv")
write.table(x = genecnv_long_df, file = file2write, quote = F, sep = "\t", row.names = F)
