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
exp_wide_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/preprocess/filter_and_transform_DIA_protein_data/20210205.v1/RCC_PDX.DIA_Protein.Log2.Filtered.20210205.v1.tsv")
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")
## input the DEGs
data_status_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/0.sample_info/v3.20210116/RCC_PDX_Samples.20210115.v2.xlsx")

# identify genes to plot -------------------------------------------------
gene_plot_df <- data.frame(gene_symbol = c("HGF", "VEGFA", "EFNA5", "EFNB2", "TGFA", "NCAM1", "NCAM1", "IL34"))

# for (i in 1:nrow(gene_plot_df)) {
for (i in c(1)) {
  gene_plot <- gene_plot_df$gene_symbol[i]
  # make plot data ----------------------------------------------------------
  plotdata_wide_df <- exp_wide_df %>%
    filter(PG.Genes %in% gene_plot)
  if (nrow(plotdata_wide_df) == 0) {
    plotdata_wide_df <- exp_wide_df %>%
      filter(grepl(x = PG.Genes, pattern = gene_plot))
  }
  if (nrow(plotdata_wide_df) == 0) {
   next() 
  }
  plotdata_df <- melt(data = plotdata_wide_df)
  plotdata_df <- plotdata_df %>%
    rename(SampleID = variable)
  plotdata_df <- merge(x = plotdata_df, y = sampleinfo_df, by.x = c("SampleID"),  by.y = c("Sample ID"), all.x = T)
  ## set cap values
  plotdata_df <- plotdata_df %>%
    mutate(ModelID = str_split_fixed(string = SampleID, pattern = "_", n = 3)[,1]) %>%
    filter(ModelID %in% c("RESL5", "RESL10")) %>%
    mutate(x_plot = value) %>%
    filter(!is.na(x_plot))
  plotdata_df$y_plot <- mapvalues(x = plotdata_df$Treatment, from = c("Con", "Cabo", "Cabo+ Sap", "Sap"), to = c("CT", "Cabo", "Cabo+Sap", "Sap"))
  plotdata_df$y_plot <- factor(x = plotdata_df$y_plot, levels = rev(c("CT", "Cabo", "Sap", "Cabo+Sap")))
  ## make colors
  colors_model <- RColorBrewer::brewer.pal(n = 6, name = "Set1")[c(1,2)]
  names(colors_model) <- unique(plotdata_df$ModelID)
  # plot --------------------------------------------------------------------
  p <- ggplot(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot, color = ModelID))
  p <- p + geom_point(alpha = 0.7, size = 3)
  p <- p + scale_color_manual(values = colors_model)
  p <- p + facet_grid(rows = vars(Treatment_length))
  p <- p + theme_classic(base_size = 12)
  p <- p + xlab("bulk protein level (log2 intensity)")
  p <- p + theme(axis.text.y = element_text(size = 12), axis.title.y = element_blank())
  p <- p + scale_x_reverse() + scale_y_discrete(position = "right")
  p <- p + theme(legend.position = "left")
  p <- p + theme(axis.text.x = element_text(size = 12), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "first")))
  p <- p + ggtitle(label = paste0(gene_plot, " bulk protein abundance (Human)"))
  file2write <- paste0(dir_out, gene_plot, ".bysample", ".png")
  png(file2write, width = 700, height = 350, res = 150)
  print(p)
  dev.off()
}
