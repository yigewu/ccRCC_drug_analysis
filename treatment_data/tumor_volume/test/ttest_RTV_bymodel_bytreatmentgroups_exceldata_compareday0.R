# Yige Wu @WashU March 2022

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
packages = c(
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "readxl"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg_name_tmp)
  }
  library(package = pkg_name_tmp, character.only = T)
}


# input data --------------------------------------------------
## input the RTV data
rtv_df <- read_xlsx(path = "./Resources/Treatment_Lists/Treatment_Summary/Relative_Tumor_Volume_tests.xlsx", sheet = "RTV")

# process -----------------------------------------------------------------
rtv_df <- rtv_df %>%
  mutate(id_model_batch = paste0(Model, "_B", Batch)) %>%
  filter(Treatment_days > 0)
colnames_rtv <- colnames(rtv_df); colnames_rtv <- colnames_rtv[grepl(pattern = "RTV", x = colnames_rtv)]
treatment_groups <- unique(rtv_df$Treatment_group)
test_res_df <- NULL

for (id_model_batch_tmp in unique(rtv_df$id_model_batch)) {
  rtv_tmp_df <- rtv_df %>%
    filter(id_model_batch == id_model_batch_tmp)
  for (day_tmp in unique(rtv_tmp_df$Treatment_days)) {
    rtv_tmp2_df <- rtv_tmp_df %>%
      filter(Treatment_days == day_tmp)
    wilcox_pvalue_vec <- NULL
    wilcox_w_vec <- NULL
    ttest_pvalue_vec <- NULL
    ttest_t_vec <- NULL
    for (treatment_group in treatment_groups) {
      rtv_group1 <- rtv_tmp2_df[rtv_tmp2_df$Treatment_group == treatment_group, colnames_rtv]; rtv_group1 <- rtv_group1[!is.na(rtv_group1)]
      if(length(rtv_group1) > 1) {
        rtv_group2 = rep(100, length(rtv_group1))
        wilcox_res <- wilcox.test(x = rtv_group1, y = rtv_group2)
        wilcox_pvalue_vec <- c(wilcox_pvalue_vec, wilcox_res$p.value)
        wilcox_w_vec <- c(wilcox_w_vec, wilcox_res$statistic)
        
        ttest_res <- t.test(x = rtv_group1, y = rtv_group2)
        ttest_pvalue_vec <- c(ttest_pvalue_vec, ttest_res$p.value)
        ttest_t_vec <- c(ttest_t_vec, ttest_res$statistic)
      } else {
        wilcox_pvalue_vec <- c(wilcox_pvalue_vec, NA)
        wilcox_w_vec <- c(wilcox_w_vec, NA)
        ttest_pvalue_vec <- c(ttest_pvalue_vec, NA)
        ttest_t_vec <- c(ttest_t_vec, NA)
      }
    }
    test_res_tmp_df <- data.frame(treatment_groups, ttest_pvalue_vec, ttest_t_vec, wilcox_pvalue_vec, wilcox_w_vec)
    colnames(test_res_tmp_df) <- c("Treatment_group", "Pvalue_ttest", "T_ttest", "Pvalue_wilcox", "W_wilcox")
    test_res_tmp_df$FDR_ttest <- p.adjust(test_res_tmp_df$Pvalue_ttest, method = "fdr")
    test_res_tmp_df$id_model_batch <- id_model_batch_tmp
    test_res_tmp_df$Treatment_days <- day_tmp
    test_res_df <- rbind(test_res_df, test_res_tmp_df)
  }
}

# write output ------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_drug_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
file2write <- paste0(dir_out, "Treatment_result_test.", run_id, ".tsv")
write.table(x = test_res_df, file = file2write, sep = "\t", row.names = F, quote = F)

