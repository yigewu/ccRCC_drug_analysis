# Yige Wu @WashU Apr 2022
## plot cell type on integration UMAP

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2",
  "ggrastr"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
source("./ccRCC_drug_analysis//functions.R")
## set run id
version_tmp <- 3
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
plot_data_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_byres_humancells_8sample_integration_on_katmai/20220905.v1/HumanCells.8sample.Metadata.ByResolution.20220905.v1.tsv")

# preprocess --------------------------------------------------------------
plot_data_df <- plot_data_df %>%
  mutate(model = str_split_fixed(string = orig.ident, pattern = "\\-", n = 3)[,1]) %>%
  mutate(model = gsub(x = model, pattern = '[A-Z]', replacement = "")) %>%
  mutate(model = paste0("RESL", model)) %>%
  mutate(treatment = str_split_fixed(string = orig.ident, pattern = "\\-", n = 3)[,3]) %>%
  mutate(treatment = gsub(x = treatment, pattern = '[0-9]', replacement = ""))
table(plot_data_df$treatment)
plot_data_df$treatment[plot_data_df$treatment == "CT"] <- "control" 
plot_data_df$treatment[plot_data_df$treatment == "Sap"] <- "sapanisertib" 
plot_data_df$treatment[plot_data_df$treatment == "Cabo"] <- "cabozantinib" 
plot_data_df$treatment[plot_data_df$treatment == "Cabo_Sap"] <- "cabozantinib+sapanisertib" 
plot_data_df$treatment <- factor(x = plot_data_df$treatment, levels = c("control", "cabozantinib", "sapanisertib", "cabozantinib+sapanisertib"))
## adjust transparency by cell number
cellnumber_df <- plot_data_df %>%
  group_by(orig.ident) %>%
  summarise(cellnumber = n())
cellnumber_df$alpha = log(min(cellnumber_df$cellnumber))/log(cellnumber_df$cellnumber)
plot_data_df$alpha <- mapvalues(x = plot_data_df$orig.ident, from = cellnumber_df$orig.ident, to = as.vector(cellnumber_df$alpha))
plot_data_df$alpha <- as.numeric(plot_data_df$alpha)
plot_data_df$size <- 0.1*plot_data_df$alpha
# make plots --------------------------------------------------------------
for (res_tmp in c(0.5)) {
# for (res_tmp in c(0.1, 0.2, 0.3, 0.4, 0.5, 1, 2)) {
  plot_data_df[, "cluster_id"] <- plot_data_df[, paste0("integrated_snn_res.", res_tmp)]
  plot_data_df$clustername_plot <- paste0("MC", (plot_data_df$cluster_id+1))
  clusterids <- paste0("MC", as.numeric(sort(unique(factor(plot_data_df$cluster_id)))))
  colors_bycluster <- Polychrome::palette36.colors(n = (length(clusterids)+1))[-2]
  names(colors_bycluster) <- clusterids
  p <- ggplot()
  # p <- p + geom_point_rast(data = plot_data_df, 
  #                          mapping = aes(x = UMAP_1, y = UMAP_2, color = clustername_plot),
  #                          alpha = 1, size = 0.1, shape = 16)
  # p <- p + geom_point_rast(data = plot_data_df, 
  #                          mapping = aes(x = UMAP_1, y = UMAP_2, color = clustername_plot, alpha = alpha),
  #                          size = 0.1, shape = 16)
  p <- p + geom_point_rast(data = plot_data_df, 
                           mapping = aes(x = UMAP_1, y = UMAP_2, color = clustername_plot, size = alpha),
                           alpha = 0.7, shape = 16) +
    scale_size_continuous(range = c(0.05, 0.7))
  p <- p + facet_grid(rows = vars(model), cols = vars(treatment))
  p <- p + scale_color_manual(values = colors_bycluster)
  p <- p + guides(colour = guide_legend(override.aes = list(size=4), title = NULL, label.theme = element_text(size = 18), nrow = 3))
  p <- p + theme_void()
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_blank())
  # axis.line = element_line(colour = "black"))
  p <- p + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())
  p <- p + theme(axis.text.y=element_blank(),
                 axis.ticks.y=element_blank())
  p <- p + theme(legend.position="none", aspect.ratio=1)
  ## save as png
  file2write <- paste0(dir_out, "resolution", res_tmp, ".png")
  png(filename = file2write, width = 1000, height = 1100, res = 150)
  print(p)
  dev.off()
  file2write <- paste0(dir_out, "resolution", res_tmp, ".pdf")
  pdf(file2write, width = 10, height = 6, useDingbats = F)
  print(p)
  dev.off()
}

