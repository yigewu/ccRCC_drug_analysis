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
barcode2cluster_df <- fread(input = "./Resources/Analysis_Results/snrna_processing/map_barcode/map_barcode2idmanualcluster_humancells_integration_withanchor/20210212.v1/barcode2idmanualcluster_umapdata.20210212.v1.tsv", data.table = F)

# make plot data ----------------------------------------------------------
plot_data_df <- barcode2cluster_df %>%
  mutate(orig.ident = Id_Sample) %>%
  group_by(orig.ident, Id_Manual_Cluster) %>%
  summarise(count_bysample_bycluster = n())
count_bysample_df <- barcode2cluster_df %>%
  mutate(orig.ident = Id_Sample) %>%
  group_by(orig.ident) %>%
  summarise(count_bysample = n())  
plot_data_df <- merge(x = plot_data_df, y = count_bysample_df, by = c("orig.ident"), all.x = T)
plot_data_df <- plot_data_df %>%
  mutate(frac_bysample_bycluster = (count_bysample_bycluster/count_bysample)) %>%
  mutate(id_model = str_split_fixed(string = orig.ident, pattern = "-", n = 3)[,1]) %>%
  mutate(treatment_group = str_split_fixed(string = orig.ident, pattern = "-", n = 3)[,3]) %>%
  mutate(treatment_group = gsub(pattern = "2", replacement = "", x = treatment_group)) %>%
  mutate(treatment_group = factor(x = treatment_group, levels = c("CT", "Sap", "Cabo", "Cabo_Sap"))) %>%
  mutate(cluster = factor(Id_Manual_Cluster))
## make color for each cluster
colors_cluster <- Polychrome::dark.colors(n = length(unique(plot_data_df$Id_Manual_Cluster)))
names(colors_cluster) <- unique(plot_data_df$Id_Manual_Cluster)

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
p <- p + facet_grid(cols = vars(id_model), scales = "free")
p <- p + theme_classic()
p <- p + guides(colour = guide_legend(ncol = 2))
p <- p + ylim(c(0, 0.1))
p
file2write <- paste0(dir_out, "CT-Sap-Combo.", "limit10", ".png")
png(filename = file2write, width = 1000, height = 600, res = 150)
print(p)
dev.off()

# make plot for CT-Sap-Combo ---------------------------------------------------------
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

p <- ggplot(data = subset(plot_data_df, treatment_group != "Sap"), mapping = aes(x = treatment_group, y = frac_bysample_bycluster, group = cluster, color = cluster))
p <- p + geom_line()
p <- p + geom_point()
p <- p + scale_color_manual(values = colors_cluster)
p <- p + facet_grid(cols = vars(id_model), scales = "free")
p <- p + theme_classic()
p <- p + guides(colour = guide_legend(ncol = 2))
p <- p + ylim(c(0, 0.1))
p
file2write <- paste0(dir_out, "CT-Cabo-Combo.", "limit10",  ".png")
png(filename = file2write, width = 1000, height = 600, res = 150)
print(p)
dev.off()

