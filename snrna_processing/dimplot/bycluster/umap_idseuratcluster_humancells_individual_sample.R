# Yige Wu @WashU Feb 2021
## plot dimplot with cluster name

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
## input cell type and umap data
barcode2cluster_df <- fread(input = "./Resources/Analysis_Results/snrna_processing/fetch_data/fetch_data_individual_sample_humancells_on_katmai/20210224.v1/HumanCells.Individual_sample.umap_data.20210224.v1.tsv", data.table = F)

# make umap by sample  ---------------------------------------------------------------
for (sampleid_tmp in unique(barcode2cluster_df$orig.ident)) {
  ## make plot data
  plotdata_df <- barcode2cluster_df %>%
    filter(orig.ident == sampleid_tmp) %>%
    mutate(id_cluster = factor(seurat_clusters))
  ## make color for each cluster
  colors_cluster <- Polychrome::dark.colors(n = length(unique(plotdata_df$seurat_clusters)))
  names(colors_cluster) <- unique(plotdata_df$seurat_clusters)
  ## make plot
  p <- ggplot()
  p <- p + geom_point(data = plotdata_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = id_cluster), 
                      size = 0.2, alpha = 0.8)
  p <- p + scale_color_manual(values = colors_cluster)
  p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p <- p + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())
  p <- p + theme(axis.text.y=element_blank(),
                 axis.ticks.y=element_blank())
  p
  ## write output
  file2write <- paste0(dir_out, "umap.", sampleid_tmp, ".", "png")
  png(filename = file2write, width = 1200, height = 1000, res = 150)
  print(p)
  dev.off()
}



