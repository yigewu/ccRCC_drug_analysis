# Yige Wu @ WashU 2023 Jun

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
packages = c(
  "plyr",
  "stringr",
  "reshape2",
  "data.table",
  "dplyr"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

# input dependencies ------------------------------------------------------
de_df <- fread("./Resources/Analysis_Results/expression/rna/test/ttest_paired_diff_TFactivity_1month_combo_vs_single/20230621.v4/TFactivity.Ttest.Paired.1month.Combo_vs_single.20230621.v4.tsv", data.table = F)
## input ERK downstream
genes2filter_df = readxl::read_xlsx(path = "./Resources/Knowledge/41580_2020_255_MOESM1_ESM.xlsx", skip = 1)
target2gene_df = OmnipathR::import_transcriptional_interactions()

# filter ------------------------------------------------------------------
genes_select_df = exp_test_df %>%
  filter(group2 == "Control") %>%
  filter(Name %in% genes2filter_df$SUBSTRATE_GENE) %>%
  filter(diff_estimate < 0) %>% ## Combo < Control
  filter(!(Name %in% c("FOXO3"))) %>%
  arrange(pvalue)

genes_filter1_df = exp_test_df %>%
  filter(group2 == "Treated.Cabo") %>%
  filter(Name %in% genes2filter_df$SUBSTRATE_GENE) %>%
  filter(diff_estimate < 0)

genes_filter2_df = exp_test_df %>%
  filter(group2 == "Treated.Sap") %>%
  filter(Name %in% genes2filter_df$SUBSTRATE_GENE) %>%
  filter(diff_estimate < 0)

genes_select_df = genes_select_df %>%
  filter(Name %in% genes_filter1_df$Name) %>%
  filter(Name %in% genes_filter2_df$Name)

target2gene_filtered_df = target2gene_df %>%
  filter(source_genesymbol %in% genes_select_df$Name) %>%
  filter(target_genesymbol %in% c("CCND3", "MYBL2", "FOXM1", "CDK2", "BCL2L12", "BCL2L13"))

# save output -------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_drug_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
