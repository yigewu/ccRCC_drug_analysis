# Yige Wu @ WashU 2021 Nov

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
packages = c(
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table"
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

# input dependencies ------------------------------------------------------
result_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/rna/test/ttest_paired_diff_test_1month_combo_vs_single/20220322.v1/mRNA.Ttest.Paired.1month.Combo_vs_single.20220322.v1.tsv")
## input druggable genes
druggable_df1 <- fread(data.table = F, input = "./Resources/Knowledge/Gene_Lists/approved_target_ids_all.csv")
druggable_df2 <- readxl::read_excel(path = "../ccRCC_snRNA/Resources/Knowledge/Gene_Lists/Targetable_Genes.20200924.xlsx")

# pre-processing ---------------------------------------------
genes_targetable <- unique(c(druggable_df1$`Gene Name`, druggable_df2$genesymbol))

plotdata_df <- result_df %>%
  filter(group2 == "Control") %>%
  # filter(pvalue < 0.05 & diff_estimate >= 0.1) %>%
  filter(pvalue < 0.05 & diff_estimate > 0) %>%
  filter(Name %in% genes_targetable) %>%
  mutate(x_plot = diff_estimate) %>%
  arrange(diff_estimate)
