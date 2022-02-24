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
exp_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/preprocess/filter_and_transform_DIA_phosphorylation_data/20210205.v1/RCC_PDX.DIA_Phosphopeptide.Log2.20210205.v1.tsv")
## input the sample info
meta_data_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")
## input phosphosite location
ptm_names_byprotein_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/preprocess/locate_phosphosites/20210204.v1/Phosphorylation_Sites.20210204.v1.tsv")

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
exp_process_df <- merge(x = exp_process_df, y = ptm_names_byprotein_df %>%
                          select(PG.ProteinName, EG.ModifiedSequence, PTM_Name),
                        by = c("PG.ProteinName", "EG.ModifiedSequence"),
                        all.x = T)

# run by loop -------------------------------------------------------------
meta_data_control_df <- meta_data_df %>%
  filter(Treatment == "Con") %>%
  mutate(ModelID = str_split_fixed(string = `Sample ID`, pattern = "_", n = 3)[,1]) %>%
  arrange(ModelID)
meta_data_treated_df <- meta_data_df %>%
  filter(Treatment != "Con") %>%
  mutate(ModelID = str_split_fixed(string = `Sample ID`, pattern = "_", n = 3)[,1]) %>%
  arrange(ModelID)

# model_id_tmp <- "RESL10"
# treatment_tmp <- "Cabo"
# treatment.month_tmp <- "1 month"
foldchange_all_df <- exp_process_df[, c("PG.Gene", "PG.ProteinName", "PG.ProteinDescriptions", "PG.ProteinAccession", "PTM_Name")]

for (model_id_tmp in unique(meta_data_control_df$ModelID)) {
  meta_data_treated_tmp_df <- meta_data_treated_df %>%
    filter(ModelID == model_id_tmp)
  for (treatment_tmp in unique(meta_data_treated_tmp_df$Treatment)) {
    meta_data_treated_tmp2_df <- meta_data_treated_tmp_df %>%
      filter(Treatment == treatment_tmp)
    for (treatment.month_tmp in unique(meta_data_treated_tmp2_df$Treatment_length)) {
      analysis_ids_tmp <- meta_data_control_df$`Sample ID`[meta_data_control_df$ModelID == model_id_tmp & meta_data_control_df$Treatment_length == treatment.month_tmp]
      exp_control_vec <- exp_process_df[, analysis_ids_tmp]
      
      analysis_ids_tmp <- meta_data_treated_df$`Sample ID`[meta_data_treated_df$ModelID == model_id_tmp & meta_data_treated_df$Treatment == treatment_tmp & meta_data_treated_df$Treatment_length == treatment.month_tmp]
      exp_treated_vec <- exp_process_df[, analysis_ids_tmp]
      
      foldchange_vec <- (exp_treated_vec - exp_control_vec)
      foldchange_df <- data.frame(foldchange_vec); colnames(foldchange_df) <- paste0(model_id_tmp, "_", treatment_tmp, "_", gsub(x = treatment.month_tmp, pattern = " ", replacement = ""))
      foldchange_all_df <- cbind(foldchange_all_df, foldchange_df)
    }
  }
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Log2Phosphorylation_diff_treated_vs_control.", run_id, ".tsv")
write.table(x = foldchange_all_df, file = file2write, quote = F, sep = "\t", row.names = F)
