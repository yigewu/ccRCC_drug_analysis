# Yige Wu @ WashU 2023 Jun

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("OmnipathR")
packages = c(
  "plyr",
  "stringr",
  "reshape2",
  "data.table",
  "dplyr",
  "ggplot2"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input dependencies ------------------------------------------------------
## input the protein data
protein_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/preprocess/normalize_DIA_protein_data/20210111.v1/RCC_PDX.DIA_Protein.Log2.QuantileNormalized.20210111.v1.tsv")
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")

# set parameters ----------------------------------------------------------
genes_filter <- c("ICAM1")
genes_filter <- c("CCND3", "MYBL2", "FOXM1", "CDK2", "BCL2L12","BCL2L13")
genes_filter = c("JUNB", "ELK1", "SMAD4", "BRF1", "POU5F1", "FOSL2", "JUND", "FOSL1", "ERG",
  "FOS", "SRF", "EGR1", "ELK4", "JUN", "ETV1", "MYC", "MITF", "STAT3")

colnames_id <- colnames(protein_df)[!(colnames(protein_df) %in% sampleinfo_df$`Sample ID`)]
genes_filter = genes_filter[genes_filter %in% protein_df$PG.Genes]; genes_filter

# plot --------------------------------------------------------------------
for (gene_plot in genes_filter) {
  ## make plot data
  plot_data_wide_df <- protein_df %>%
    filter(PG.Genes == gene_plot)
  plot_data_long_df <- melt(data = plot_data_wide_df, id.vars = colnames_id)
  plot_data_long_df$Treatment_length <- mapvalues(x = plot_data_long_df$variable, from = sampleinfo_df$`Sample ID`, to = as.vector(sampleinfo_df$Treatment_length))
  plot_data_long_df$Treatment <- mapvalues(x = plot_data_long_df$variable, from = sampleinfo_df$`Sample ID`, to = as.vector(sampleinfo_df$Treatment))
  plot_data_long_df <- as.data.frame(plot_data_long_df)
  plot_data_long_df$y_plot <- (plot_data_long_df$value - min(plot_data_long_df$value, na.rm = T))/(max(plot_data_long_df$value, na.rm = T) - min(plot_data_long_df$value, na.rm = T))
  plot_data_long_df <- plot_data_long_df %>%
    dplyr::mutate(Model_id = str_split_fixed(string = variable, pattern = "_", n = 3)[,1]) %>%
    dplyr::filter(Treatment %in% c("Con", "Cabo+ Sap")) %>%
    dplyr::filter(Treatment_length == "1 month")
  
  # plot_data_long_df$Treatment <- factor(x = plot_data_long_df$Treatment, levels = c("Con", "Sap", "Cabo", "Cabo+ Sap"))
  plot_data_long_df$Treatment <- factor(x = plot_data_long_df$Treatment, levels = c("Con", "Cabo", "Sap", "Cabo+ Sap"))
  plot_data_long_df$Model_id <- factor(x =   plot_data_long_df$Model_id, levels = c("RESL5", "RESL10", "RESL11", "RESL3", "RESL12", "RESL4"))
  
  test_df = plot_data_long_df %>%
    filter(Model_id != "RESL11")
  # test_df = test_df %>%
  #   filter(Model_id != "RESL5")
  print(t.test(test_df$y_plot[test_df$Treatment == "Cabo+ Sap"], test_df$y_plot[test_df$Treatment == "Con"], paired = T))
  
  ## plot
  p <- ggline(data = plot_data_long_df, 
                 x = "Treatment", y = "y_plot", color = "Model_id",
                 # add = c("mean", "dotplot"),
                 position = position_dodge())
  # p <- p + ylim(c(-0.05, 1.05))
  p <- p + labs(y = paste0(gene_plot, " protein level\n(normalized)"))
  p <- p + theme_classic(base_size = 15)
  p <- p + theme(#axis.text.x = element_blank(),
                 # axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(),
                 legend.position = "right")
  p
  ## write output
  file2write <- paste0(dir_out, gene_plot, ".pdf")
  pdf(file2write, width = 4, height = 3, useDingbats = F)
  print(p)
  dev.off()
  file2write <- paste0(dir_out, gene_plot, ".png")
  png(file2write, width = 1500, height = 800, res = 150)
  print(p)
  dev.off()
}
