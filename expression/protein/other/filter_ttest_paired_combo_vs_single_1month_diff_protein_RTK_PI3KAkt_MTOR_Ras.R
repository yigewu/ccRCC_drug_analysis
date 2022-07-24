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

# input pathway member ----------------------------------------------------
pathwaynames_ordered <- c("RTK_ligand", "PI3K_Akt", "TSC_mTOR", "Ras_MAPK")
pathway2members_df <- NULL
for (pathwayname_plot in pathwaynames_ordered) {
  pathway2members_tmp_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Gene_Lists/Pathway_score_members.011222.xlsx", sheet = pathwayname_plot)
  pathway2members_tmp_df$Pathway_name <- pathwayname_plot
  pathway2members_df <- rbind(pathway2members_df, pathway2members_tmp_df[, c("Name", "Gene_symbol", "Is_phospho", "Site_phospho", "Type_of_regulation", "Pathway_name")])
}
pathway2members_df <- pathway2members_df %>%
  mutate(ID = paste0(Gene_symbol, "_", ifelse(Is_phospho == "Yes", Site_phospho, "Protein")))

# filter for same direction in both comparisons-----------------------------------------------
result_filtered_df <-  result_merged_df %>%
  filter(group2 == "Cabo") %>%
  mutate(diff_direction = ifelse(diff_estimate > 0, "up", "down")) %>%
  mutate(gene_species = ifelse(grepl(pattern = "HUMAN", PG.ProteinNames), 
                               ifelse(grepl(pattern = "MOUSE", PG.ProteinNames), "Ambiguous", "Human"), "Mouse")) %>%
  filter(PG.Genes %in% pathway2members_df$Gene_symbol) %>%
  arrange(pvalue)

result_filtered_df <-  result_merged_df %>%
  filter(group2 == "Sap") %>%
  mutate(diff_direction = ifelse(diff_estimate > 0, "up", "down")) %>%
  mutate(gene_species = ifelse(grepl(pattern = "HUMAN", PG.ProteinNames), 
                               ifelse(grepl(pattern = "MOUSE", PG.ProteinNames), "Ambiguous", "Human"), "Mouse")) %>%
  filter(PG.Genes %in% pathway2members_df$Gene_symbol) %>%
  arrange(pvalue)


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
