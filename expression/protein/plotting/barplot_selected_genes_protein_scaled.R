# Yige Wu @ WashU 2021 Mar

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
source("./ccRCC_drug_analysis/plotting.R")
## set run id
version_tmp <- 1
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
# genes_filter <- c("MT-CO1", "MT-CO3")
# genes_filter <- c("TMEM63A", "Aig1;AIG1", "C11orf54", "SLC27A3")
genes_filter <- c("ICAM1")
colnames_id <- colnames(protein_df)[!(colnames(protein_df) %in% sampleinfo_df$`Sample ID`)]

# plot --------------------------------------------------------------------
gene_plot <- genes_filter[1]
for (gene_plot in genes_filter) {
  ## make plot data
  plot_data_wide_df <- protein_df %>%
    filter(PG.Genes == gene_plot)
  plot_data_long_df <- melt(data = plot_data_wide_df, id.vars = colnames_id)
  plot_data_long_df$Treatment_length <- mapvalues(x = plot_data_long_df$variable, from = sampleinfo_df$`Sample ID`, to = as.vector(sampleinfo_df$Treatment_length))
  plot_data_long_df$Treatment <- mapvalues(x = plot_data_long_df$variable, from = sampleinfo_df$`Sample ID`, to = as.vector(sampleinfo_df$Treatment))
  plot_data_long_df <- as.data.frame(plot_data_long_df)
  plot_data_long_df <- plot_data_long_df %>%
    dplyr::mutate(Model_id = str_split_fixed(string = variable, pattern = "_", n = 3)[,1])
  # plot_data_long_df$Treatment <- factor(x = plot_data_long_df$Treatment, levels = c("Con", "Sap", "Cabo", "Cabo+ Sap"))
  plot_data_long_df$Treatment <- factor(x = plot_data_long_df$Treatment, levels = c("Con", "Cabo", "Sap", "Cabo+ Sap"))
  plot_data_long_df$y_plot <- scale(x = plot_data_long_df$value)
  plot_data_long_df$Model_id <- factor(x =   plot_data_long_df$Model_id, levels = c("RESL5", "RESL4", "RESL10", "RESL3", "RESL11", "RESL12", "RESL6", "RESL8", "RESL9"))
  ## plot
  p <- ggplot()
  p <- p + geom_bar(data = plot_data_long_df, mapping = aes(x = Treatment, y = y_plot, fill = Model_id), stat = "identity")
  p <- p + geom_hline(yintercept = 0, linetype = 2)
  p <- p + facet_grid(cols = vars(Model_id), rows = vars(Treatment_length), scales = "free")
  p <- p + theme_classic()
  p <- p + theme(axis.text.x = element_text(angle = 15))
  # p
  ## write output
  file2write <- paste0(dir_out, gene_plot, ".png")
  png(file2write, width = 1500, height = 800, res = 150)
  print(p)
  dev.off()
}
