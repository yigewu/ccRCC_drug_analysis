# Yige Wu @WashU May 2020
## plot dimplot with cluster name

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
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
## input cell type and umap data
barcode2cluster_df <- fread(input = "./Resources/Analysis_Results/snrna_processing/fetch_data/fetch_data_humancells_8sample_integration_withanchor_on_katmai/20210209.v1/HumanCells_8sample.umap_data.20210209.v1.tsv", data.table = F)

# make plot data ----------------------------------------------------------
plot_data_df <- barcode2cluster_df %>%
  group_by(orig.ident, seurat_clusters) %>%
  summarise(count_bysample_bycluster = n())
count_bysample_df <- barcode2cluster_df %>%
  group_by(orig.ident) %>%
  summarise(count_bysample = n())  
plot_data_df <- merge(x = plot_data_df, y = count_bysample_df, by = c("orig.ident"), all.x = T)
plot_data_df <- plot_data_df %>%
  mutate(frac_bysample_bycluster = (count_bysample_bycluster/count_bysample)) %>%
  mutate(id_model = str_split_fixed(string = orig.ident, pattern = "-", n = 3)[,1]) %>%
  mutate(treatment_group = str_split_fixed(string = orig.ident, pattern = "-", n = 3)[,3]) %>%
  mutate(treatment_group = gsub(pattern = "2", replacement = "", x = treatment_group)) %>%
  mutate(treatment_group = factor(x = treatment_group, levels = c("CT", "Sap", "Cabo", "Cabo_Sap"))) %>%
  mutate(cluster = paste0("MC", (seurat_clusters+1))) %>%
  mutate(cluster = factor(cluster))
## make color for each cluster
clusterids <- paste0("MC", as.numeric(sort(unique(factor(plot_data_df$seurat_clusters)))))
colors_cluster <- Polychrome::palette36.colors(n = (length(clusterids)+1))[-2]
names(colors_cluster) <- clusterids

# make plot for CT-Sap-Combo ---------------------------------------------------------
p <- ggplot(data = subset(plot_data_df, treatment_group != "Cabo"), mapping = aes(x = treatment_group, y = frac_bysample_bycluster, group = cluster, color = cluster))
p <- p + geom_line()
p <- p + geom_point()
p <- p + scale_color_manual(values = colors_cluster)
p <- p + facet_grid(cols = vars(id_model), scales = "free")
p <- p + theme_classic()
p <- p + guides(colour = guide_legend(ncol = 2))

p
file2write <- paste0(dir_out, "CT-Sap-Combo", ".png")
png(filename = file2write, width = 1000, height = 600, res = 150)
print(p)
dev.off()

p <- ggplot(data = subset(plot_data_df, treatment_group != "Cabo"), mapping = aes(x = treatment_group, y = frac_bysample_bycluster, group = cluster, color = cluster))
p <- p + geom_line()
p <- p + geom_point()
p <- p + scale_color_manual(values = colors_cluster)
p <- p + facet_grid(rows = vars(id_model), cols = vars(cluster), scales = "free")
p <- p + theme_classic()
p <- p + guides(colour = guide_legend(ncol = 2))

file2write <- paste0(dir_out, "CT-Sap-Combo.split", ".png")
png(filename = file2write, width = 2000, height = 1000, res = 150)
print(p)
dev.off()

# make plot for CT-Cab-Combo ---------------------------------------------------------
p <- ggplot(data = subset(plot_data_df, treatment_group != "Sap"), mapping = aes(x = treatment_group, y = frac_bysample_bycluster, group = cluster, color = cluster))
p <- p + geom_line()
p <- p + geom_point()
p <- p + scale_color_manual(values = colors_cluster)
p <- p + facet_grid(cols = vars(id_model), scales = "free")
p <- p + theme_classic()
p <- p + guides(colour = guide_legend(ncol = 2))

p
file2write <- paste0(dir_out, "CT-Cabo-Combo", ".png")
png(filename = file2write, width = 1000, height = 600, res = 150)
print(p)
dev.off()

