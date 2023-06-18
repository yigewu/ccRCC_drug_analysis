# Yige Wu @WashU May 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the average expression
exp_wide_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/average_expression/unite_expression/unite_sct_data_humancells_mousecells_bycelltypeshorter_bysample/20210505.v1/SCT.data.AverageExpression.ByCellTypeShorter.BySample.20210505.v1.tsv")

# identify genes to plot -------------------------------------------------
x_cap <- 20
genes_plot <- c("PXDN", "PDXP", "GLUD2", "SUSD2")
gene_plot <- "PXDN"
genes_plot <- c("C3", "C1R")

## make colors
colors_cellgroup <- RColorBrewer::brewer.pal(n = 5, name = "Set1")
names(colors_cellgroup) <- c("Tumor.cells", "Endothelial.cells", "Macrophages", "Fibroblasts", "Myofibroblasts")

for (gene_plot in genes_plot) {
  # make plot data ----------------------------------------------------------
  plotdata_wide_df <- exp_wide_df %>%
    filter(genesymbol_human %in% gene_plot)
  plotdata_df <- melt(data = plotdata_wide_df)
  summary(plotdata_df$value)
  plotdata_df <- plotdata_df %>%
    mutate(cell_group = str_split_fixed(string = variable, pattern = "2_", n = 2)[,2]) %>%
    mutate(sample = str_split_fixed(string = variable, pattern = "2_", n = 2)[,1]) %>%
    mutate(treatment = str_split_fixed(string = sample, pattern = "\\.", n = 3)[,3]) %>%
    mutate(model = str_split_fixed(string = sample, pattern = "\\.", n = 3)[,1]) %>%
    mutate(y_plot = paste0(model, "-", treatment)) %>%
    mutate(x_plot = ifelse(value > x_cap, x_cap, value)) %>%
    arrange(model, factor(x = treatment, levels = c("CT", "Cabo", "Sap", "Cabo_Sap")))
  
  plotdata_df$y_plot <- factor(x = plotdata_df$y_plot, levels = rev(unique(plotdata_df$y_plot)))

  # plot --------------------------------------------------------------------
  p <- ggplot()
  p <- p + geom_point(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot, color = cell_group), alpha = 0.7, size = 3)
  p <- p + scale_color_manual(values = colors_cellgroup)
  p <- p + theme_classic(base_size = 12)
  p <- p + xlab("Normalized expression")
  p <- p + theme(axis.text.y = element_text(size = 12), axis.title.y = element_blank())
  p <- p + theme(axis.text.x = element_text(size = 12))
  p <- p + scale_x_reverse() + scale_y_discrete(position = "right")
  p <- p + theme(legend.position = "left")
  file2write <- paste0(dir_out, gene_plot, ".bysample", ".png")
  png(file2write, width = 1000, height = 500, res = 150)
  print(p)
  dev.off()
}
