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
result_merged_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/test/ttest_paired_diff_protein_1month_combo_vs_single/20220316.v1/Ttest.Paired.1month.Combo_vs_single.20220316.v1.tsv")

# see the distribution of the protein difference --------------------------
ggplot() +
  geom_histogram(data = subset(x = result_merged_df, group2 == "Cabo"), mapping = aes(x = diff_estimate))
result_merged_df %>%
  filter(group2 == "Cabo") %>%
  nrow()
## [1] 7520

result_merged_df %>%
  filter(group2 == "Cabo") %>%
  filter(pvalue < 0.05) %>%
  nrow()
## [1] 353

result_merged_df %>%
  filter(group2 == "Cabo") %>%
  filter(pvalue < 0.05) %>%
  filter(abs(diff_estimate) >= 0.1) %>%
  nrow()
## [1] 352

# result_merged_df %>%
#   filter(group2 == "Cabo") %>%
#   filter(pvalue < 0.05) %>%
#   filter(abs(diff_estimate) >= 0.1) %>%
#   View()

result_merged_df %>%
  filter(group2 == "Sap") %>%
  nrow()
## [1] 7829

result_merged_df %>%
  filter(group2 == "Sap") %>%
  filter(pvalue < 0.05) %>%
  nrow()
## [1] 430

result_merged_df %>%
  filter(group2 == "Sap") %>%
  filter(pvalue < 0.05) %>%
  filter(abs(diff_estimate) >= 0.1) %>%
  nrow()
## [1] 427

# result_merged_df %>%
#   filter(group2 == "Cabo") %>%
#   filter(pvalue < 0.05) %>%
#   filter(abs(diff_estimate) >= 0.1) %>%
#   View()

# filter by p-value and diff------------------------------------------------------------------
cutoff_diff <- 0.1
cutoff_pvalue <- 0.05
result_filtered_df <- result_merged_df %>%
  filter(pvalue < cutoff_pvalue) %>%
  filter(abs(diff_estimate) >= cutoff_diff) %>%
  arrange(pvalue)

# filter for same direction in both comparisons-----------------------------------------------
result_filtered_df %>%
  filter(group2 == "Cabo") %>%
  mutate(diff_direction = ifelse(diff_estimate > 0, "up", "down")) %>%
  mutate(gene_species = ifelse(grepl(pattern = "HUMAN", PG.ProteinNames), 
                               ifelse(grepl(pattern = "MOUSE", PG.ProteinNames), "Ambiguous", "Human"), "Mouse")) %>%
  select(PG.Genes, diff_direction, gene_species) %>%
  unique() %>%
  select(diff_direction, gene_species) %>%
  table()
# gene_species
# diff_direction Ambiguous Human Mouse
# down        11   146    14
# up          19    59   103

result_filtered_df %>%
  filter(group2 == "Sap") %>%
  mutate(diff_direction = ifelse(diff_estimate > 0, "up", "down")) %>%
  mutate(gene_species = ifelse(grepl(pattern = "HUMAN", PG.ProteinNames), 
                               ifelse(grepl(pattern = "MOUSE", PG.ProteinNames), "Ambiguous", "Human"), "Mouse")) %>%
  select(PG.Genes, diff_direction, gene_species) %>%
  unique() %>%
  select(diff_direction, gene_species) %>%
  table()
# gene_species
# diff_direction Ambiguous Human Mouse
# down        34   154   165
# up           7    42    23

intersect(subset(x = result_filtered_df, group2 == "Cabo")$PG.Genes, subset(x = result_filtered_df, group2 == "Sap")$PG.Genes)
intersect(subset(x = result_filtered_df, group2 == "Cabo" & diff_estimate > 0)$PG.Genes, subset(x = result_filtered_df, group2 == "Sap" & diff_estimate > 0)$PG.Genes)
intersect(subset(x = result_filtered_df, group2 == "Cabo" & diff_estimate < 0)$PG.Genes, subset(x = result_filtered_df, group2 == "Sap" & diff_estimate < 0)$PG.Genes)

# save output -------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_drug_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
file2write <- paste0(dir_out, "Ttest.Paired.1month.Combo_vs_single.P", cutoff_pvalue, ".Diff", cutoff_diff, ".", run_id, ".tsv")
write.table(x = result_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
