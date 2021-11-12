# Yige Wu @WashU May 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
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
gene2cellgroup_plot_df <- data.frame(gene_symbol = c("FLT1", "FLT1", "KDR", "MET", "EPHA3", "EGFR", "EGFR", "PDGFRB", "IGF1R", "IGF1R", "IGF1R", "CSF1R", "FGFR1"),
                                     cell_group = c("Endothelial.cells", "Myofibroblasts", "Endothelial.cells", "Tumor.cells", "Fibroblasts", "Tumor.cells", "Fibroblasts", "Myofibroblasts", "Tumor.cells", "Fibroblasts", "Myofibroblasts", "Macrophages", "Fibroblasts"))

gene_plot <- genes_plot[1]
cellgroup_plot <- "Tumor.cells"
for (i in 1:nrow(gene2cellgroup_plot_df)) {
  # for (i in c(1)) {
    
  gene_plot <- gene2cellgroup_plot_df$gene_symbol[i]
  cellgroup_plot <- gene2cellgroup_plot_df$cell_group[i]
  # make plot data ----------------------------------------------------------
  plotdata_wide_df <- exp_wide_df %>%
    filter(genesymbol_human %in% gene_plot)
  plotdata_df <- melt(data = plotdata_wide_df)
  ## set cap values
  plotdata_df <- plotdata_df %>%
    mutate(cell_group = str_split_fixed(string = variable, pattern = "2_", n = 2)[,2]) %>%
    filter(cell_group == cellgroup_plot) %>%
    mutate(sample = str_split_fixed(string = variable, pattern = "2_", n = 2)[,1]) %>%
    mutate(treatment = str_split_fixed(string = sample, pattern = "\\.", n = 3)[,3]) %>%
    mutate(model = str_split_fixed(string = sample, pattern = "\\.", n = 3)[,1]) %>%
    mutate(x_plot = value) %>%
    mutate(treatment_length = "2-month")
  plotdata_df$y_plot <- factor(x = plotdata_df$treatment, levels = rev(c("CT", "Cabo", "Sap", "Cabo_Sap")))
  ## make colors
  colors_model <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[c(1,2)]
  names(colors_model) <- unique(plotdata_df$model)
  # plot --------------------------------------------------------------------
  p <- ggplot()
  p <- p + geom_point(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot, color = model), alpha = 0.7, size = 3)
  p <- p + scale_color_manual(values = colors_model)
  p <- p + facet_grid(rows = vars(treatment_length))
  p <- p + theme_classic(base_size = 12)
  p <- p + xlab("Normalized expression")
  p <- p + theme(axis.text.y = element_text(size = 12), axis.title.y = element_blank())
  p <- p + theme(axis.text.x = element_text(size = 12))
  p <- p + ggtitle(label = paste0(gene_plot, " sn RNA expression", ifelse(cellgroup_plot == "Tumor.cells", " (Human)", " (Mouse)")),
                   subtitle = cellgroup_plot)
  file2write <- paste0(dir_out, gene_plot, ".", cellgroup_plot, ".bysample", ".png")
  png(file2write, width = 750, height = 275, res = 150)
  print(p)
  dev.off()
}
