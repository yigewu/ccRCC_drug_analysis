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
gi_group_sup_df <- fread(data.table = F, input = "./Resources/Analysis_Results/treatment_data/extract_cellline_growth_inhibtion_bydrugs_083021/20210903.v1/RCC_Cellline.Growth_Inhibition.20210903.v1.tsv")

# write by drug -----------------------------------------------------------
drug_tmp <- "MK571"
cellline_tmp <- "Hk2"

for (drug_tmp in unique(gi_group_sup_df$Drug[!is.na(gi_group_sup_df$Drug)])) {
  dir_out_tmp <- paste0(dir_out, drug_tmp, "/")
  dir.create(dir_out_tmp)
  for (cellline_tmp in unique(gi_group_sup_df$cellline)) {
    gi_tmp_df <- gi_group_sup_df %>%
      filter(cellline == cellline_tmp) %>%
      filter(Drug == drug_tmp) %>%
      mutate(repname = paste0("rep", repid))
    gi_ctl_tmp_df <- gi_group_sup_df %>%
      filter(cellline == cellline_tmp) %>%
      filter(drug_token %in% c("negcntrl", "pstvcntrl")) %>%
      mutate(repname = paste0("rep", ifelse(conc_token != "Media", repid, ceiling(repid/2)))) %>%
      mutate(conc_value = paste0(drug_token, "_", conc_token))
    gi_merged_tmp_df <- rbind(gi_tmp_df, gi_ctl_tmp_df)
    gi_merged_tmp2_df <- gi_merged_tmp_df %>%
      group_by(conc_value, repname) %>%
      summarise(growth_inhibition = mean(growth_inhibition))
    gi_wide_tmp_df <- dcast(data = gi_merged_tmp2_df, formula = conc_value ~ repname, value.var = "growth_inhibition")
    
    file2write <- paste0(dir_out_tmp, drug_tmp, ".", cellline_tmp, ".Growth_inhibition.tsv")
    write.table(x = gi_wide_tmp_df, file = file2write, sep = "\t", row.names = F, quote = F)
  }
}
