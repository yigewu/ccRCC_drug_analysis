# Yige Wu @WashU Jan 2022

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
  "ggplot2",
  "ggrastr",
  "ggrepel"
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
## input DEGs
enricher_top_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/pathway/unite_ora_msigdb_H_CP_treated_vs_control_human_proteins/20220216.v1/ora_msigdb_H_CP.treated_vs_control.human.proteins.20220216.v1.tsv")
## input test result
result_merged_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/test/ttest_paired_diff_protein_1month_treated_vs_control/20220323.v1/Ttest.Paired.1month.Treated_vs_control.20220323.v1.tsv")

# preprocess -------------------------------------------------------------
genes_highlight_tmp <- c("TK1", "CDK1", "KPNA2", "MCM2", "PRDX4", "PCNA", "TBRG4", "MCM3", "RFC3", "MCM7",
                         #"TIMM8B", "GOT2", "HCCS", "UQCRQ", "SLC25A20", "NDUFB1", "COX7C", "TOMM22", "SLC25A5", "UQCRB", "UQCR10", "COX4I1", "COX17", "MRPL11", "DLD")#,
                          "SCG2", "LRP1", "TGFBI", "THBS1", "COL7A1", "TFPI2", "TAGLN",
                          "CDK4", "TNPO2", "RPS6KA5", "PRMT5", "NUP98")
                          #"ICAM1", "LRP1", "FTH1", "ICAM1", "ALDH1L1", "CDK4", "BLVRB", "ADK", "MAT1A")
cap_value <- 2
## make plot data
plot_data_df <- result_merged_df %>%
  filter(group1 == "Cabo+ Sap") %>%
  mutate(gene_species = ifelse(grepl(pattern = "HUMAN", PG.ProteinNames), 
                               ifelse(grepl(pattern = "MOUSE", PG.ProteinNames), "Ambiguous", "Human"), "Mouse")) %>%
  filter(gene_species == "Human") %>%
  mutate(diff_direction = ifelse(diff_estimate > 0, "up", "down")) %>%
  mutate(text_gene = PG.Genes) %>%
  mutate(y_plot = -log10(pvalue)) %>%
  # mutate(x_plot = diff_estimate) %>%
  mutate(x_plot = ifelse(diff_estimate > cap_value, cap_value,
                         ifelse(diff_estimate < (-cap_value), (-cap_value), diff_estimate))) %>%
  mutate(group = ifelse(pvalue < 0.05 & abs(diff_estimate) > 0.1, ifelse(diff_estimate > 0, "up", "down"), "insignificant")) %>%
  unique()

plot_data_df <- plot_data_df %>%
  mutate(label_plot = ifelse(text_gene %in% genes_highlight_tmp & pvalue < 0.05, text_gene, NA))

# plot ------------------------------------------------------------------
fontsize_plot = 12
p <- ggplot(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot, color = group, label = label_plot))
p <- p + geom_point_rast(alpha = 0.6, size = 1, shape = 16)
p <- p + geom_text_repel(data = subset(plot_data_df, x_plot > 0),
                         force = 5, xlim = c(0.8, 2.75),
                         segment.size = 0.2, segment.alpha = 0.5,
                         max.overlaps = Inf)
p <- p + geom_text_repel(data = subset(plot_data_df, x_plot < 0),
                         force = 5, xlim = c(-2.55, -1),
                         segment.size = 0.2, segment.alpha = 0.5,
                         max.overlaps = Inf)
p <- p + xlim(-2.5, 2.5)
p <- p + scale_color_manual(values = c("up" = "#E41A1C", "down" = "#377EB8", "insignificant" = "grey80"))
p <- p + xlab("Difference in log2(Intensity) (combination-treated - control)") + ylab("-Log10(P-value)")
p <- p + ggtitle(label = paste0("Protein changes btw. cab.+sap. vs. vehicle treatment"))
p <- p + theme_bw() + theme(legend.position = "none", title = element_text(size = 10), 
                            axis.text = element_text(size = fontsize_plot), axis.title = element_text(size = fontsize_plot))
# p <- p + ggtitle(label = paste0(sampleid1, " vs ", sampleid2, "\nMouse Endothelial Cells snRNA Expression"))

# save output -------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_drug_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
file2write <- paste0(dir_out, "volcano", ".pdf")
pdf(file2write, width = 5, height = 4, useDingbats = F)
print(p)
dev.off()
# file2write <- paste0(dir_out, treatment_tmp, ".png")
# png(filename = file2write, width = 800, height = 600, res = 150)
# print(p)
# dev.off()
