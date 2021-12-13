# Yige Wu @ WashU 2021 Nov
## need to overlap with protein changes

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
## input fold changes
fc_rna_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/rna/get_foldchange/get_foldchange_treated_vs_baseline/20211109.v1/RNA_foldchange_treated_vs_baseline.20211109.v1.tsv")
fc_protein_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/get_foldchange_treated_vs_control_byprotein/20211109.v1/Protein_diff_treated_vs_control.20211109.v1.tsv")
## input detailed sample meta data
meta_data_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/0.sample_info/v3.20210116/RCC_PDX_Samples.20210115.v2.xlsx")

# run by loop -------------------------------------------------------------
meta_data_baseline_df <- meta_data_df %>%
  filter(ShortTag == "Baseline") %>%
  filter(DataType == "RNA-Seq") %>%
  arrange(ModelID)
meta_data_treated_df <- meta_data_df %>%
  filter(ShortTag %in% c(treatment_tmp, "Control")) %>%
  filter(DataType == "RNA-Seq") %>%
  arrange(ModelID)

# process sensitivity-related genes -----------------------------------------------------------
fc_protein_filtered_df <- fc_protein_df %>%
  filter(RESL5_Cabo_1month < 0 & RESL5_Cabo_1month < RESL5_Sap_1month)
fc_rna_filtered_df <- fc_rna_df %>%
  filter(!is.infinite(RESL5_Treated.Cabo_1month) & RESL5_Treated.Cabo_1month < 1 & 
           (RESL5_Treated.Cabo_1month/RESL5_Control_1month) < 1 & 
           (RESL5_Treated.Cabo_1month < RESL5_Treated.Sap_1month))
nrow(fc_rna_filtered_df)
fc_rna_filtered_df <- fc_rna_filtered_df %>%
  filter(genesymbol %in% fc_protein_filtered_df$PG.Gene)
nrow(fc_rna_filtered_df)
fc_rna_filtered_df$gene_category <- "Cabo sensitive"
fn_rna_filtered_merged_df <- fc_rna_filtered_df

# process resistance-related genes -----------------------------------------------------------
fc_protein_filtered_df <- fc_protein_df %>%
  filter(RESL10_Cabo_1month > 0 & RESL4_Cabo_1month > 0 & RESL10_Cabo_1month > RESL10_Sap_1month & RESL4_Cabo_1month > RESL4_Sap_1month)
fc_rna_filtered_df <- fc_rna_df %>%
  filter(RESL10_Treated.Cabo_1month > 0 & !is.infinite(RESL10_Treated.Cabo_1month) & 
           (RESL4_Treated.Cabo_1month >= 1 | RESL10_Treated.Cabo_1month >= 1) & 
           (RESL10_Treated.Cabo_1month/RESL10_Control_1month) >= 1 & (RESL10_Treated.Cabo_1month/RESL10_Treated.Sap_1month) >= 1 & 
           (RESL4_Treated.Cabo_1month/RESL4_Control_1month) >= 1 & (RESL4_Treated.Cabo_1month/RESL4_Treated.Sap_1month) >= 1)
nrow(fc_rna_filtered_df)
fc_rna_filtered_df <- fc_rna_filtered_df %>%
  filter(genesymbol %in% fc_protein_filtered_df$PG.Gene)
nrow(fc_rna_filtered_df)
fc_rna_filtered_df$gene_category <- "Cabo resistant"
fn_rna_filtered_merged_df <- rbind(fn_rna_filtered_merged_df, fc_rna_filtered_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Cabo_related_genes.", run_id, ".tsv")
write.table(x = fn_rna_filtered_merged_df,  file = file2write, quote = F, sep = "\t", row.names = F)

# make background gene list ----------------------------------------------
fc_rna_df <- data.frame(fc_rna_df)
background_genes_df <- fc_rna_df %>%
  filter(!is.infinite(RESL5_Treated.Cabo_1month)) %>%
  filter(genesymbol %in% fc_protein_df$PG.Gene[!is.na(fc_protein_df$RESL5_Cabo_1month) & !is.na(fc_protein_df$RESL5_Sap_1month)]) %>%
  dplyr::select(genesymbol)
file2write <- paste0(dir_out, "Cabo_sensitive_background_genes.", run_id, ".tsv")
write.table(x = background_genes_df, file = file2write, quote = F, sep = "\t", row.names = F)

# make background gene list ----------------------------------------------
background_genes_df <- fc_rna_df %>%
  filter(!is.infinite(RESL10_Treated.Cabo_1month) & !is.infinite(RESL4_Treated.Cabo_1month)) %>%
  filter(genesymbol %in% fc_protein_df$PG.Gene[!is.na(fc_protein_df$RESL10_Cabo_1month) & 
                                                 !is.na(fc_protein_df$RESL4_Cabo_1month) & 
                                                 !is.na(fc_protein_df$RESL10_Sap_1month) & 
                                                 !is.na(fc_protein_df$RESL4_Sap_1month)]) %>%
  dplyr::select(genesymbol)
file2write <- paste0(dir_out, "Cabo_resistant_background_genes.", run_id, ".tsv")
write.table(x = background_genes_df, file = file2write, quote = F, sep = "\t", row.names = F)

