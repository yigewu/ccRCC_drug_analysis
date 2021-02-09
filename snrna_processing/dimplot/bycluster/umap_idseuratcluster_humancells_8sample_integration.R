# Yige Wu @WashU May 2020
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
barcode2cluster_df <- fread(input = "./Resources/Analysis_Results/snrna_processing/fetch_data/fetch_data_humancells_8sample_integration_withanchor_on_katmai/20210209.v1/HumanCells_8sample.umap_data.20210209.v1.tsv", data.table = F)

# make umap with all samples & cluster together ---------------------------------------------------------------
## make plot data
plotdata_df <- barcode2cluster_df %>%
  mutate(id_seurat_cluster = factor(seurat_clusters))
## make color for each cluster
colors_cluster <- Polychrome::dark.colors(n = length(unique(plotdata_df$id_seurat_cluster)))
names(colors_cluster) <- unique(plotdata_df$id_seurat_cluster)
## make plot
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = id_seurat_cluster), 
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
file2write <- paste0(dir_out, "umap.", "allclusters.", "allsamples.",  ".png")
png(filename = file2write, width = 1200, height = 1000, res = 150)
print(p)
dev.off()

# make umap by RESL with all clusters together ---------------------------------------------------------------
## make plot data
plotdata_df <- barcode2cluster_df %>%
  mutate(id_seurat_cluster = factor(seurat_clusters)) %>%
  mutate(id_model = str_split_fixed(string = orig.ident, pattern = "-", n = 3)[,1])
## make plot
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = id_seurat_cluster), 
                    size = 0.2, alpha = 0.8)
p <- p + scale_color_manual(values = colors_cluster)
# p <- p + facet_grid(.~orig.ident, rows = vars(id_model))
p <- p + facet_grid(cols = vars(id_model))
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p
## write output
file2write <- paste0(dir_out, "umap.", "allclusters.", "byRESL.",  ".png")
png(filename = file2write, width = 2200, height = 1000, res = 150)
print(p)
dev.off()

# make umap by sample with all clusters together ---------------------------------------------------------------
## make plot data
plotdata_df <- barcode2cluster_df %>%
  mutate(id_seurat_cluster = factor(seurat_clusters)) %>%
  mutate(id_model = str_split_fixed(string = orig.ident, pattern = "-", n = 3)[,1]) %>%
  mutate(treatment_group = str_split_fixed(string = orig.ident, pattern = "-", n = 3)[,3]) %>%
  mutate(treatment_group = gsub(pattern = "2", replacement = "", x = treatment_group))
plotdata_df$treatment_group <- factor(x = plotdata_df$treatment_group, levels = c("CT", "Sap", "Cabo", "Cabo_Sap"))
## make plot
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = id_seurat_cluster), 
                    size = 0.2, alpha = 0.8)
p <- p + scale_color_manual(values = colors_cluster)
# p <- p + facet_grid(.~orig.ident, rows = vars(id_model))
p <- p + facet_grid(rows = vars(id_model), cols = vars(treatment_group))
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p
## write output
file2write <- paste0(dir_out, "umap.", "allclusters.", "bysample.",  ".png")
png(filename = file2write, width = 2200, height = 1000, res = 150)
print(p)
dev.off()

# make umap by cluster with all samples together ---------------------------------------------------------------
## set output directory
dir_out1 <- paste0(dir_out, "ByCluster_AllSamples", "/")
dir.create(dir_out1)
for (cluster_tmp in unique(barcode2cluster_df$seurat_clusters)) {
  ## make plot data
  plotdata_df <- barcode2cluster_df %>%
    mutate(id_seurat_cluster = ifelse(seurat_clusters == cluster_tmp, cluster_tmp, "others"))
  ## make color for each cluster
  colors_cluster_tmp <- c(colors_cluster[as.character(cluster_tmp)], "grey80")
  names(colors_cluster_tmp) <- c(cluster_tmp, "others")
  ## make plot
  p <- ggplot()
  p <- p + geom_point(data = plotdata_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = id_seurat_cluster), 
                      size = 0.2, alpha = 0.8)
  p <- p + scale_color_manual(values = colors_cluster_tmp)
  p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p <- p + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())
  p <- p + theme(axis.text.y=element_blank(),
                 axis.ticks.y=element_blank())
  p
  ## write output
  file2write <- paste0(dir_out1, "umap.", "clusters", cluster_tmp, ".", "allsamples.",  ".png")
  png(filename = file2write, width = 1200, height = 1000, res = 150)
  print(p)
  dev.off()
}



