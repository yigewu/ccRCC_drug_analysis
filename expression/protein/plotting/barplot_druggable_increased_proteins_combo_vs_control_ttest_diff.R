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
  "data.table",
  "ggplot2"
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
result_filtered_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/test/ttest_paired_diff_protein_1month_combo_vs_single/20220317.v1/Ttest.Paired.1month.Combo_vs_single.20220317.v1.tsv")
## input druggable genes
druggable_df1 <- fread(data.table = F, input = "./Resources/Knowledge/Gene_Lists/approved_target_ids_all.csv")
druggable_df2 <- readxl::read_excel(path = "../ccRCC_snRNA/Resources/Knowledge/Gene_Lists/Targetable_Genes.20200924.xlsx")
## input the sample info
meta_data_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")

# prepare plotting parameters ---------------------------------------------
cap_value <- 2
# cap_value <- 7

genes_targetable <- unique(c(druggable_df1$`Gene Name`, druggable_df2$genesymbol))
# genes_targetable[!(genes_targetable %in% result_filtered_df$PG.Genes[result_filtered_df$group2 == "Con"])]
druggable_df2$genesymbol[!(druggable_df2$genesymbol %in% result_filtered_df$PG.Genes[result_filtered_df$group2 == "Con"])]

plotdata_df <- result_filtered_df %>%
  filter(group2 == "Con") %>%
  # filter(pvalue < 0.05 & diff_estimate >= 0.1) %>%
  filter(pvalue < 0.05 & diff_estimate > 0.1) %>%
  mutate(gene_species = ifelse(grepl(pattern = "HUMAN", PG.ProteinNames), 
                               ifelse(grepl(pattern = "MOUSE", PG.ProteinNames), "Ambiguous", "Human"), "Mouse")) %>%
  
  filter(gene_species == "Human") %>%
  filter(PG.Genes %in% genes_targetable) %>%
  mutate(x_plot = diff_estimate) %>%
  arrange(diff_estimate)
plotdata_df$y_plot <- factor(x = plotdata_df$PG.Genes, levels = plotdata_df$PG.Genes)

## plot
p <- ggplot()
# p <- p + geom_vline(xintercept = 0, linetype = 2)
p <- p + geom_bar(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot), stat = "identity", fill = "#E41A1C")
# p <- p + scale_fill_manual(values = colors_treatment)
p <- p + theme_classic()
p <- p + xlab("Increase in log2(Intensity)\n(combination-treated - control)") 
p <- p + theme(axis.title.y = element_blank())
p

# save output -------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_drug_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
file2write <- paste0(dir_out, "allmodels", "_", "1month", ".pdf")
pdf(file2write, width = 4, height = 2.5, useDingbats = F)
print(p)
dev.off()


