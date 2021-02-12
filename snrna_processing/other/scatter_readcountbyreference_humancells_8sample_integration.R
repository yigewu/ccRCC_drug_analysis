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
## input readcount info
readcount_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/filtering/make_seuratfiltered_srat_obj_katmai/20210127.v1/CellRanger_Filtered.Barcode_Metrics.tsv")

# specify parameters ------------------------------------------------------
## make color for each cluster
colors_cluster <- Polychrome::dark.colors(n = length(unique(barcode2cluster_df$seurat_clusters)))
names(colors_cluster) <- unique(barcode2cluster_df$seurat_clusters)


# preprocess --------------------------------------------------------------
plotdata_df <- merge(x = barcode2cluster_df %>%
                       select(-V1) %>%
                       mutate(barcode_raw = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1]), 
                     y = readcount_df, by = c("orig.ident", "barcode_raw"), all.x = T)
plotdata_df <- plotdata_df %>%
  mutate(id_seurat_cluster = factor(seurat_clusters)) %>%
  mutate(id_model = str_split_fixed(string = orig.ident, pattern = "-", n = 3)[,1]) %>%
  mutate(treatment_group = str_split_fixed(string = orig.ident, pattern = "-", n = 3)[,3]) %>%
  mutate(treatment_group = gsub(pattern = "2", replacement = "", x = treatment_group))
plotdata_df$treatment_group <- factor(x = plotdata_df$treatment_group, levels = c("CT", "Sap", "Cabo", "Cabo_Sap"))

# plot for all clusters all samples ---------------------------------------------------------------
## make plot
p <- ggplot()
p <- p + geom_abline(slope = 1, linetype = 2)
p <- p + geom_hline(yintercept = 4000, linetype = 2)

p <- p + geom_point(data = plotdata_df, mapping = aes(x = mm10, y = GRCh38, color = id_seurat_cluster), 
                    size = 0.2, alpha = 0.8)
p <- p + scale_color_manual(values = colors_cluster)
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + theme_classic()
p
## write output
file2write <- paste0(dir_out, "allclusters.", "allsamples.",  ".png")
png(filename = file2write, width = 1200, height = 1000, res = 150)
print(p)
dev.off()

# plot for all clusters all samples -same x, y axis scale--------------------------------------------------------------
## make plot
p <- ggplot()
p <- p + geom_abline(slope = 1, linetype = 2)
p <- p + geom_point(data = plotdata_df, mapping = aes(x = mm10, y = GRCh38, color = id_seurat_cluster), 
                    size = 0.2, alpha = 0.8)
p <- p + xlim(c(0, 4000)) + ylim(c(0, 4000))
p <- p + scale_color_manual(values = colors_cluster)
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + theme_classic()
p
## write output
file2write <- paste0(dir_out, "allclusters.", "allsamples.samescale.",  ".png")
png(filename = file2write, width = 1200, height = 1000, res = 150)
print(p)
dev.off()

# plot for by clusters all samples ---------------------------------------------------------------
## set output directory
dir_out1 <- paste0(dir_out, "ByCluster_AllSamples", "/")
dir.create(dir_out1)
for (cluster_tmp in unique(barcode2cluster_df$seurat_clusters)) {
  ## make plot data
  plotdata_df_tmp <- plotdata_df %>%
    mutate(id_seurat_cluster = ifelse(seurat_clusters == cluster_tmp, cluster_tmp, "others"))
  ## make color for each cluster
  colors_cluster_tmp <- c(colors_cluster[as.character(cluster_tmp)], "grey80")
  names(colors_cluster_tmp) <- c(cluster_tmp, "others")
  
  ## make plot
  p <- ggplot()
  p <- p + geom_abline(slope = 1, linetype = 2)
  p <- p + geom_point(data = subset(plotdata_df_tmp, id_seurat_cluster != cluster_tmp), mapping = aes(x = mm10, y = GRCh38, color = id_seurat_cluster), 
                      size = 0.2, alpha = 0.8)
  p <- p + geom_point(data = subset(plotdata_df_tmp, id_seurat_cluster == cluster_tmp), mapping = aes(x = mm10, y = GRCh38, color = id_seurat_cluster), 
                      size = 0.2, alpha = 0.8)
  p <- p + scale_color_manual(values = colors_cluster_tmp)
  p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
  p <- p + theme_classic()
  p
  ## write output
  file2write <- paste0(dir_out1, "cluster", cluster_tmp, ".allsamples",  ".png")
  png(filename = file2write, width = 1200, height = 1000, res = 150)
  print(p)
  dev.off()
}

# plot for by clusters all samples same x, y axis scale---------------------------------------------------------------
## set output directory
dir_out1 <- paste0(dir_out, "ByCluster_AllSamples_Same_Scale", "/")
dir.create(dir_out1)
for (cluster_tmp in unique(barcode2cluster_df$seurat_clusters)) {
  ## make plot data
  plotdata_df_tmp <- plotdata_df %>%
    mutate(id_seurat_cluster = ifelse(seurat_clusters == cluster_tmp, cluster_tmp, "others"))
  ## make color for each cluster
  colors_cluster_tmp <- c(colors_cluster[as.character(cluster_tmp)], "grey80")
  names(colors_cluster_tmp) <- c(cluster_tmp, "others")
  
  ## make plot
  p <- ggplot()
  p <- p + geom_abline(slope = 1, linetype = 2)
  p <- p + geom_point(data = subset(plotdata_df_tmp, id_seurat_cluster != cluster_tmp), mapping = aes(x = mm10, y = GRCh38, color = id_seurat_cluster), 
                      size = 0.2, alpha = 0.8)
  p <- p + geom_point(data = subset(plotdata_df_tmp, id_seurat_cluster == cluster_tmp), mapping = aes(x = mm10, y = GRCh38, color = id_seurat_cluster), 
                      size = 0.2, alpha = 0.8)
  p <- p + scale_color_manual(values = colors_cluster_tmp)
  p <- p + xlim(c(0, 4000)) + ylim(c(0, 4000))
  p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
  p <- p + theme_classic()
  p
  ## write output
  file2write <- paste0(dir_out1, "cluster", cluster_tmp, ".allsamples",  ".png")
  png(filename = file2write, width = 1200, height = 1000, res = 150)
  print(p)
  dev.off()
}


