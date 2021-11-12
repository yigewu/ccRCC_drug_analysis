# Yige Wu @ WashU 2020 Jun

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
library(readxl)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
gi_group_sup_df <- fread(data.table = F, input = "./Resources/Analysis_Results/treatment_data/cell_line_treatment/extract_growth_inhibition/extract_cellline_growth_inhibtion_bycombodrugs_092321/20210924.v1/RCC_Cellline.Cabo_Sapa_TM.Combo.Growth_Inhibition.20210924.v1.tsv")
std <- function(x) sd(x)/sqrt(length(x))

# write by drug -----------------------------------------------------------
drug_tmp <- unique(gi_group_sup_df$Drug[!is.na(gi_group_sup_df$Drug)])[1]

for (drug_tmp in unique(gi_group_sup_df$Drug[!is.na(gi_group_sup_df$Drug)])) {
  gi_tmp_df <- gi_group_sup_df %>%
    filter(Drug == drug_tmp) %>%
    mutate(repname = paste0("rep", repid))
  gi_ctl_tmp_df <- gi_group_sup_df %>%
    filter(drug_token %in% c("negcntrl", "pstvcntrl")) %>%
    mutate(repname = paste0("rep", repid)) %>%
    mutate(conc_value = paste0(drug_token, "_", conc_token))
  gi_merged_tmp_df <- rbind(gi_tmp_df, gi_ctl_tmp_df)
  ## sumarize mean growth inhibition by cell line by drug
  gi_merged_tmp_df <- gi_merged_tmp_df %>%
    mutate(cell_rep = paste0(cellline, "_", repname))
  gi_wide_df <- dcast(data = gi_merged_tmp_df, formula = conc_value ~ cell_rep, value.var = "growth_inhibition")
  file2write <- paste0(dir_out, drug_tmp, ".Growthinhibition.", "ByCellLine.", "ByReplicate.tsv")
  write.table(x = gi_wide_df, file = file2write, sep = "\t", row.names = F, quote = F)
  
  ## sumarize mean growth inhibition by cell line by drug
  gi_mean_long_df <- gi_merged_tmp_df %>%
    group_by(cellline, conc_value) %>%
    summarise(growth_inhibition = mean(growth_inhibition))
  gi_mean_wide_df <- dcast(data = gi_mean_long_df, formula = conc_value ~ cellline, value.var = "growth_inhibition")
  file2write <- paste0(dir_out, drug_tmp, ".Growthinhibition.", "ByCellLine.", "MeanAcrossReplicates.tsv")
  write.table(x = gi_mean_wide_df, file = file2write, sep = "\t", row.names = F, quote = F)
  
  ## sumarize SEM for growth inhibition by cell line by dru
  gi_sem_long_df <- gi_merged_tmp_df %>%
    group_by(cellline, conc_value) %>%
    summarise(sem_growth_inhibition = std(growth_inhibition))
  gi_sem_wide_df <- dcast(data = gi_sem_long_df, formula = conc_value ~ cellline, value.var = "sem_growth_inhibition")
  file2write <- paste0(dir_out, drug_tmp, ".Growthinhibition.", "ByCellLine.", "SEMAcrossReplicates.tsv")
  write.table(x = gi_sem_wide_df, file = file2write, sep = "\t", row.names = F, quote = F)
  
}
