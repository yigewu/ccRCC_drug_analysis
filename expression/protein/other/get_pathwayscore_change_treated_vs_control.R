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
exp_process_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/calculate_pathway_score_bysample/20220110.v1/Pathway_scores_bysample.20220110.v1.tsv")
## input the sample info
meta_data_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")

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
foldchange_all_df <- data.frame(Pathway_name = exp_process_df$Pathway_name)
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
file2write <- paste0(dir_out, "Pathway_score_diff_treated_vs_control.", run_id, ".tsv")
write.table(x = foldchange_all_df, file = file2write, quote = F, sep = "\t", row.names = F)
