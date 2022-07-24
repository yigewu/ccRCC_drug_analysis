# Yige Wu @ WashU 2022 Feb

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
  "ggplot2"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

# input dependencies ------------------------------------------------------
## input the test result
result_merged_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/test/ttest_paired_diff_protein_1month_treated_vs_control/20220323.v1/Ttest.Paired.1month.Treated_vs_control.20220323.v1.tsv")

# filter by p-value and diff------------------------------------------------------------------
cutoff_diff <- 0.1
cutoff_pvalue <- 0.05
result_filtered_df <- result_merged_df %>%
  filter(pvalue < cutoff_pvalue) %>%
  filter(abs(diff_estimate) >= cutoff_diff) %>%
  mutate(gene_species = ifelse(grepl(pattern = "HUMAN", PG.ProteinNames), 
                               ifelse(grepl(pattern = "MOUSE", PG.ProteinNames), "Ambiguous", "Human"), "Mouse")) %>%
  mutate(diff_direction = ifelse(diff_estimate > 0, "up", "down")) %>%
  arrange(pvalue)

# filter for same direction in both comparisons-----------------------------------------------
result_filtered_df %>%
  filter(group1 == "Cabo+ Sap") %>%
  select(PG.Genes, diff_direction, gene_species) %>%
  unique() %>%
  select(diff_direction, gene_species) %>%
  table()
# gene_species
# diff_direction Ambiguous Human Mouse
# down        28   213    49
# up          23    97   100

result_filtered_df %>%
  filter(group1 == "Sap") %>%
  select(PG.Genes, diff_direction, gene_species) %>%
  unique() %>%
  select(diff_direction, gene_species) %>%
  table()
# gene_species
# diff_direction Ambiguous Human Mouse
# down        10   108    16
# up          37   187   213

result_filtered_df %>%
  filter(group1 == "Cabo") %>%
  select(PG.Genes, diff_direction, gene_species) %>%
  unique() %>%
  select(diff_direction, gene_species) %>%
  table()
# gene_species
# diff_direction Ambiguous Human Mouse
# down         4    43    62
# up          11    76    12

result_single_df <- result_filtered_df %>%
  filter(group1 %in% c("Cabo", "Sap")) %>%
  filter(gene_species == "Human")

intersect(subset(x = result_filtered_df, group1 == "Cabo+ Sap" & diff_estimate > 0)$PG.Genes, subset(x = result_single_df, diff_estimate > 0)$PG.Genes) %>% length()
intersect(subset(x = result_filtered_df, group1 == "Cabo+ Sap" & diff_estimate < 0)$PG.Genes, subset(x = result_single_df, diff_estimate < 0)$PG.Genes) %>% length()

# save output -------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_drug_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
file2write <- paste0(dir_out, "Ttest.Paired.1month.Treated_vs_control.P", cutoff_pvalue, ".Diff", cutoff_diff, ".", run_id, ".tsv")
write.table(x = result_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
