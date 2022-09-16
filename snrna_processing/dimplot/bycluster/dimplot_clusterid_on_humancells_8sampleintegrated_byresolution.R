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
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
plot_data_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_byres_humancells_8sample_integration_on_katmai/20220905.v1/HumanCells.8sample.Metadata.ByResolution.20220905.v1.tsv")


# make plots --------------------------------------------------------------
# for (res_tmp in c(1)) {
for (res_tmp in c(0.1, 0.2, 0.3, 0.4, 0.5, 1, 2)) {
  plot_data_df[, "cluster_id"] <- plot_data_df[, paste0("integrated_snn_res.", res_tmp)]
  plot_data_df$clustername_plot <- paste0("MC", (plot_data_df$cluster_id+1))
  clusterids <- paste0("MC", as.numeric(sort(unique(factor(plot_data_df$cluster_id)))))
  colors_bycluster <- Polychrome::palette36.colors(n = (length(clusterids)+1))[-2]
  names(colors_bycluster) <- clusterids
  p <- ggplot()
  p <- p + geom_point_rast(data = plot_data_df, 
                           mapping = aes(x = UMAP_1, y = UMAP_2, color = clustername_plot),
                           alpha = 1, size = 0.1, shape = 16)
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
  p <- p + theme(legend.position="bottom", aspect.ratio=1)
  ## save as png
  file2write <- paste0(dir_out, "resolution", res_tmp, ".png")
  png(filename = file2write, width = 1000, height = 1100, res = 150)
  print(p)
  dev.off()
  file2write <- paste0(dir_out, "resolution", res_tmp, ".pdf")
  pdf(file2write, width = 6, height = 6, useDingbats = F)
  print(p)
  dev.off()
}

