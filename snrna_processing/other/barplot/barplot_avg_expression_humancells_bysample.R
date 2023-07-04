# Yige Wu @ WashU 2023 Jun

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
packages = c(
  "plyr",
  "stringr",
  "reshape2",
  "data.table",
  "dplyr",
  "ggpubr",
  "rstatix"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

# input -------------------------------------------------------------------
exp_df = fread(data.table = F, "./Resources/Analysis_Results/snrna_processing/average_expression/avgexp_sct_data_humancells_bysample_on_katmai/20210426.v1/HumanCells.BySample.AverageExpression.20210426.v1.tsv")
genes_filter <- c("CCND3", "MYBL2", "FOXM1", "CDK2", "BCL2L12","BCL2L13")

# make data for plotting --------------------------------------------------
plot_data_wide_df = exp_df %>%
  filter(V1 %in% genes_filter)
plot_data_df = melt(plot_data_wide_df)
plot_data_df = plot_data_df %>%
  mutate(Treatment_tag = str_split_fixed(string = variable, pattern = "\\.", n = 4)[,4]) %>%
  mutate(Model_id = ifelse(grepl(x = variable, "RESL10"), "RESL10", "RESL5"))
plot_data_df$Treatment = mapvalues(x = plot_data_df$Treatment_tag,
                                   from = c("CT2", "Cabo2", "Sap2", "Cabo_Sap2"),
                                   to = c("Control", "Cabozantinib", "Sapanisertib", "Cabozantinib+Sapanisertib"))
plot_data_df$Treatment = factor(x = plot_data_df$Treatment, levels = c("Control", "Cabozantinib", "Sapanisertib", "Cabozantinib+Sapanisertib"))
plot_data_df$gene = factor(x = plot_data_df$V1, levels = genes_filter)
## make colors
color_red <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[1]
color_green <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[3]
color_yellow <- RColorBrewer::brewer.pal(n = 7, name = "Set2")[6]
color_grey <- "grey50"
colors_plot <- c("Control" = color_grey, "Cabozantinib" = color_red, "Sapanisertib" = color_green, "Cabozantinib+Sapanisertib" = color_yellow)
  
# plot --------------------------------------------------------------------
p <- ggbarplot(data = plot_data_df, 
               x = "Model_id", y = "value", fill = "Treatment",
               facet.by = "gene", nrow = 1,
               #add = c("mean"),
               position = position_dodge())
p <- p + scale_fill_manual(values = colors_plot)
p
