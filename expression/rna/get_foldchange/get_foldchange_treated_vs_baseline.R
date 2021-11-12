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
## input the protein data
rna_exp_df <- fread("./Data_Freeze/v1.dataFreeze.washU_rcc/3.geneExp/v3.20210116/datafreeze.v3.kallisto.geneExp.tsv", data.table = F)
## input detailed sample meta data
meta_data_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/0.sample_info/v3.20210116/RCC_PDX_Samples.20210115.v2.xlsx")

# run by loop -------------------------------------------------------------
meta_data_baseline_df <- meta_data_df %>%
  filter(ShortTag == "Baseline") %>%
  filter(DataType == "RNA-Seq") %>%
  arrange(ModelID)
meta_data_treated_df <- meta_data_df %>%
  filter(grepl(pattern = "Treated", x = ShortTag) | ShortTag == "Control") %>%
  filter(DataType == "RNA-Seq") %>%
  arrange(ModelID)

model_id_tmp <- "RESL10"
treatment_tmp <- "Treated.Cabo"
treatment.month_tmp <- 1
foldchange_all_df <- data.frame(genesymbol = rna_exp_df$Name)
model_ids_tmp <- unique(meta_data_baseline_df$ModelID)
model_ids_tmp <- intersect(model_ids_tmp, unique(meta_data_treated_df$ModelID))

for (model_id_tmp in model_ids_tmp) {
  analysis_ids_tmp <- meta_data_baseline_df$Analysis_ID[meta_data_baseline_df$ModelID == model_id_tmp]
  if (length(analysis_ids_tmp) > 1) {
    exp_baseline_vec <- rowMeans(log2(rna_exp_df[, analysis_ids_tmp] + 1))
  } else {
    exp_baseline_vec <- log2(rna_exp_df[, analysis_ids_tmp] + 1)
  }
  meta_data_treated_tmp_df <- meta_data_treated_df %>%
    filter(ModelID == model_id_tmp)
  for (treatment_tmp in unique(meta_data_treated_tmp_df$ShortTag)) {
    meta_data_treated_tmp2_df <- meta_data_treated_tmp_df %>%
      filter(ShortTag == treatment_tmp)
    for (treatment.month_tmp in unique(meta_data_treated_tmp2_df$Treatment.Month)) {
      analysis_ids_tmp <- meta_data_treated_df$Analysis_ID[meta_data_treated_df$ModelID == model_id_tmp & meta_data_treated_df$ShortTag == treatment_tmp & meta_data_treated_df$Treatment.Month == treatment.month_tmp]
      if (length(analysis_ids_tmp) > 1) {
        exp_treated_vec <- rowMeans(log2(rna_exp_df[, analysis_ids_tmp] + 1))
      } else {
        exp_treated_vec <- log2(rna_exp_df[, analysis_ids_tmp] + 1)
      }
      foldchange_vec <- exp_treated_vec/exp_baseline_vec
      foldchange_df <- data.frame(foldchange_vec); colnames(foldchange_df) <- paste0(model_id_tmp, "_", treatment_tmp, "_", treatment.month_tmp, "month")
      foldchange_all_df <- cbind(foldchange_all_df, foldchange_df)
    }
  }
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "RNA_foldchange_treated_vs_baseline.", run_id, ".tsv")
write.table(x = foldchange_all_df, file = file2write, quote = F, sep = "\t", row.names = F)
