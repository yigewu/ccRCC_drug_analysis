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
exp_wide_df <- fread("./Data_Freeze/v1.dataFreeze.washU_rcc/3.geneExp/v3.20210116/datafreeze.v3.kallisto.geneExp.protein_coding.tsv", data.table = F)
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")
## input the DEGs
data_status_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/0.sample_info/v3.20210116/RCC_PDX_Samples.20210115.v2.xlsx")

# identify genes to plot -------------------------------------------------
gene_plot_df <- data.frame(gene_symbol = c("FLT1", "FLT1", "KDR", "MET", "EPHA3", "EGFR", "EGFR", "PDGFRB", "IGF1R", "IGF1R", "IGF1R", "CSF1R", "FGFR1"))

# get ids in order --------------------------------------------------------
## sort data status
data_status_filtered_df <- data_status_df %>%
  filter((Analysis_ID %in% colnames(exp_wide_df))) %>%
  filter(ShortTag %in% c("Control", "Treated.Cabo", "Treated.Cabo+Sap", "Treated.Sap")) %>%
  arrange(ModelID, Treatment.Month, ShortTag)

for (i in 1:nrow(gene_plot_df)) {
# for (i in c(1)) {
  gene_plot <- gene_plot_df$gene_symbol[i]
  # make plot data ----------------------------------------------------------
  plotdata_wide_df <- exp_wide_df %>%
    filter(Name %in% gene_plot)
  plotdata_df <- melt(data = plotdata_wide_df)
  plotdata_df <- plotdata_df %>%
    rename(Analysis_ID = variable) %>%
    filter(Analysis_ID %in% data_status_filtered_df$Analysis_ID)
  plotdata_df <- merge(x = plotdata_df, y = data_status_filtered_df, by = c("Analysis_ID"), all.x = T)
  ## set cap values
  plotdata_df <- plotdata_df %>%
    filter(SampleID.AcrossDataType %in% sampleinfo_df$`Sample ID`| ModelID == "RESL3") %>%
    mutate(x_plot = value) %>%
    mutate(treatment_length = paste0(Treatment.Month, "-month"))
  plotdata_df$y_plot <- mapvalues(x = plotdata_df$ShortTag, from = c("Control", "Treated.Cabo", "Treated.Cabo+Sap", "Treated.Sap"), to = c("CT", "Cabo", "Cabo+Sap", "Sap"))
  plotdata_df$y_plot <- factor(x = plotdata_df$y_plot, levels = rev(c("CT", "Cabo", "Sap", "Cabo+Sap")))
  plotdata_df <- plotdata_df %>%
    arrange(ModelID, y_plot)
  ## make colors
  colors_model <- RColorBrewer::brewer.pal(n = 6, name = "Set1")[c(1, 3, 4, 5, 6,  2)]
  names(colors_model) <- c("RESL10", "RESL11", "RESL12", "RESL3",  "RESL4", "RESL5")
  # plot --------------------------------------------------------------------
  p <- ggplot(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot, color = ModelID))
  p <- p + geom_point(alpha = 0.7, size = 3)
  p <- p + scale_color_manual(values = colors_model)
  p <- p + facet_grid(rows = vars(treatment_length))
  p <- p + theme_classic(base_size = 12)
  p <- p + xlab("Normalized expression")
  p <- p + theme(axis.text.y = element_text(size = 12), axis.title.y = element_blank())
  p <- p + theme(axis.text.x = element_text(size = 12), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last")))
  p <- p + ggtitle(label = paste0(gene_plot, " bulk RNA expression (Human)"))
  file2write <- paste0(dir_out, gene_plot,".bysample", ".png")
  png(file2write, width = 700, height = 350, res = 150)
  print(p)
  dev.off()
}
