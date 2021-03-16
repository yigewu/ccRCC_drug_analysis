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
protein_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/preprocess/normalize_DIA_protein_data/20210111.v1/RCC_PDX.DIA_Protein.Log2.QuantileNormalized.20210111.v1.tsv")
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")

# expand by gene and filter -----------------------------------------------
## filter by NA
pro_filtered_df <- protein_df
## preprocess
pro_filtered_df[pro_filtered_df$PG.Genes == "H1-0", "PG.Genes"] <- "H1-0;H1-0"
## get repetition index
idx_rep <- sapply(pro_filtered_df$PG.Genes, function(gene_string) {
  genes_vec <- str_split(string = gene_string, pattern = ";")[[1]]
  len_rep <- length(genes_vec)
  return(len_rep)
})
genes_uniq_vec <- sapply(pro_filtered_df$PG.Genes, function(gene_string) {
  genes_vec <- str_split(string = gene_string, pattern = ";")[[1]]
  return(genes_vec)
})
proteinnames_uniq_vc <- sapply(pro_filtered_df$PG.ProteinNames, function(gene_string) {
  genes_vec <- str_split(string = gene_string, pattern = ";")[[1]]
  return(genes_vec)
})
proteinaccessions_uniq_vc <- sapply(pro_filtered_df$PG.ProteinAccessions, function(gene_string) {
  genes_vec <- str_split(string = gene_string, pattern = ";")[[1]]
  return(genes_vec)
})

## create new data frame
avg_pro_new_df <- pro_filtered_df[rep(1:nrow(pro_filtered_df), idx_rep),]
avg_pro_new_df$PG.Gene <- unlist(genes_uniq_vec)
avg_pro_new_df$PG.ProteinName <- unlist(proteinnames_uniq_vc)
avg_pro_new_df$PG.ProteinAccession <- unlist(proteinaccessions_uniq_vc)
## select columns
avg_pro_new_df <- avg_pro_new_df %>%
  mutate(log2Intensity_FC_RESL10_2m_CabvsCT = (RESL10_B1_F12467 - RESL10_B1_F12462)) %>%
  mutate(log2Intensity_FC_RESL10_2m_SapvsCT = (RESL10_B1_F12465 - RESL10_B1_F12462)) %>%
  mutate(log2Intensity_FC_RESL10_1m_CabvsCT = (RESL10_B1_F12468 - RESL10_B1_F12477)) %>%
  mutate(log2Intensity_FC_RESL10_1m_SapvsCT = (RESL10_B1_F12471 - RESL10_B1_F12477)) %>%
  mutate(log2Intensity_FC_RESL5_2m_CabvsCT = (RESL5_B2_E14537 - RESL5_B2_E14540)) %>%
  mutate(log2Intensity_FC_RESL5_2m_SapvsCT = (RESL5_B2_E14542 - RESL5_B2_E14540)) %>%
  mutate(log2Intensity_FC_RESL5_1m_CabvsCT = (RESL5_B2_E14533 - RESL5_B2_E14532)) %>%
  mutate(log2Intensity_FC_RESL5_1m_SapvsCT = (RESL5_B2_E14536 - RESL5_B2_E14532)) %>%
  select(PG.Gene, PG.ProteinName, PG.ProteinAccession, PG.ProteinNames,
         log2Intensity_FC_RESL10_2m_SapvsCT, log2Intensity_FC_RESL10_1m_SapvsCT, 
         log2Intensity_FC_RESL10_2m_CabvsCT, log2Intensity_FC_RESL10_1m_CabvsCT, 
         log2Intensity_FC_RESL5_2m_SapvsCT, log2Intensity_FC_RESL5_1m_SapvsCT,
         log2Intensity_FC_RESL5_2m_CabvsCT, log2Intensity_FC_RESL5_1m_CabvsCT)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Treated_vs_CT.", "Protein.", run_id, ".tsv")
write.table(x = avg_pro_new_df, file = file2write, quote = F, sep = "\t", row.names = F)
