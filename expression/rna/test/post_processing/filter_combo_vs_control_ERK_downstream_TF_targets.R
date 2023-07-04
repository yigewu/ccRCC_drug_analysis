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
  "clusterProfiler",
  "dplyr"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

# input dependencies ------------------------------------------------------
# de_df <- fread("./Resources/Analysis_Results/expression/rna/test/ttest_paired_diff_rna_1month_combo_vs_single/20230621.v1/mRNA.Ttest.Paired.1month.Combo_vs_single.20230621.v1.tsv", data.table = F)
de_df <- fread("./Resources/Analysis_Results/expression/rna/test/ttest_paired_diff_rna_1month_combo_vs_single/20230621.v2/mRNA.Ttest.Paired.1month.Combo_vs_single.20230621.v2.tsv", data.table = F)
## input pathway-gene
p2gene1 <- read.gmt("./Resources/Knowledge/Databases/MSigDB/h.all.v2023.1.Hs.symbols.gmt")
p2gene2 <- read.gmt("./Resources/Knowledge/Databases/MSigDB/c2.cp.v2023.1.Hs.symbols.gmt")
p2gene <- rbind(p2gene1, p2gene2)
## input ERK downstream
ERK_downstream_df = readxl::read_xlsx(path = "./Resources/Knowledge/41580_2020_255_MOESM1_ESM.xlsx", skip = 1)
## input TF-target
tf_target_df = OmnipathR::import_transcriptional_interactions()

# get genes to filter -----------------------------------------------------
## get genes that are in pathways related to cell proliferation, survival and growth
p2gene_filtered = p2gene %>%
  filter(grepl(x = term, pattern = "CELL_CYCLE|_ERK_|APOPTO"))
unique(p2gene_filtered$term)
genes_pathway = unique(p2gene_filtered$gene)

## filter pathway genes to those with known TF upstream
tf_target_filtered_df = tf_target_df %>%
  filter(source_genesymbol %in% ERK_downstream_df$SUBSTRATE_GENE)
targets_filtered = tf_target_filtered_df$target_genesymbol[tf_target_filtered_df$is_inhibition != 1]
genes_check = genes_pathway[genes_pathway %in% targets_filtered]
length(genes_check)
c("CCND1", "CCND2", "CCND3", "BCL2") %in% genes_check
genes_check = unique(c(genes_check, "JUNB", "FOSL2", "JUND", "JUN", "FOSL1", "EGR1", "FOS", "MYC"))
genes_check = unique(c("JUNB", "ELK1", "SMAD4", "BRF1", "POU5F1", "FOSL2", "JUND", "FOSL1", "ERG",
                       "FOS", "SRF", "EGR1", "ELK4", "JUN", "ETV1", "MYC", "MITF", "STAT3"))

# check their changes -----------------------------------------------------
genes_check_de_df = de_df %>%
  filter(Name %in% genes_check) %>%
  filter(group2 == "Control") %>%
  filter(diff_estimate < 0) %>%
  # filter(pvalue < 0.1) %>%
  arrange(pvalue)

genes_filter1_df = de_df %>%
  filter(Name %in% genes_check) %>%
  filter(group2 == "Treated.Cabo") %>%
  filter(diff_estimate < 0)

genes_filter2_df = de_df %>%
  filter(Name %in% genes_check) %>%
  filter(group2 == "Treated.Sap") %>%
  filter(diff_estimate < 0)

genes_filtered_df = genes_check_de_df %>%
  filter(Name %in% genes_filter1_df$Name) %>%
  filter(Name %in% genes_filter2_df$Name)

genes_filtered_df$Name
