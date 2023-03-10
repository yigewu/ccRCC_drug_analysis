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
barcode2cluster_df <- fread(input = "./Resources/Analysis_Results/snrna_processing/map_barcode/map_barcode2idmanualcluster_RESL_8sample_tumorcells_integration_withanchor/20200508.v1/barcode2idmanualcluster_umapdata.20200508.v1.tsv", data.table = F)
## set id integration and clustering
id_integration <- "RESL.Tumor_cells.integration.withanchor.20200507.v1.PC40.Res0.5"

# make umap with all samples & cluster together ---------------------------------------------------------------
## make plot data
plotdata_df <- barcode2cluster_df %>%
  mutate(Id_Manual_Cluster = factor(Id_Manual_Cluster))
## make color for each cluster
colors_cluster <- Polychrome::dark.colors(n = length(unique(plotdata_df$Id_Manual_Cluster)))
names(colors_cluster) <- unique(plotdata_df$Id_Manual_Cluster)
## make plot
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = Id_Manual_Cluster), 
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
file2write <- paste0(dir_out, "umap.", "allclusters.", "allsamples.", id_integration, ".png")
png(filename = file2write, width = 1200, height = 1000, res = 150)
print(p)
dev.off()

# make umap by RESL with all clusters together ---------------------------------------------------------------
## make plot data
plotdata_df <- barcode2cluster_df %>%
  mutate(Id_Manual_Cluster = factor(Id_Manual_Cluster)) %>%
  mutate(id_model = str_split_fixed(string = Id_Sample, pattern = "-", n = 3)[,1])
## make plot
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = Id_Manual_Cluster), 
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
file2write <- paste0(dir_out, "umap.", "allclusters.", "byRESL.", id_integration, ".png")
png(filename = file2write, width = 2200, height = 1000, res = 150)
print(p)
dev.off()

# make umap by RESL by Sample---------------------------------------------------------------
## make plot data
plotdata_df <- barcode2cluster_df %>%
  mutate(Id_Manual_Cluster = factor(Id_Manual_Cluster)) %>%
  mutate(id_model = str_split_fixed(string = Id_Sample, pattern = "-", n = 3)[,1]) %>%
  mutate(Treatment = str_split_fixed(string = Id_Sample, pattern = "-", n = 3)[,3]) %>%
  mutate(Treatment = gsub(x = Treatment, pattern = '[0-9]', replacement = ""))
plotdata_df$Treatment <- factor(x = plotdata_df$Treatment, levels = c("CT", "Cabo", "Sap", "Cabo_Sap"))
plotdata_df$id_model <- factor(x = plotdata_df$id_model, levels = c("RESL5E", "RESL10F"))

## make plot
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = Id_Manual_Cluster), 
                    size = 0.2, alpha = 0.8)
p <- p + scale_color_manual(values = colors_cluster)
# p <- p + facet_grid(.~orig.ident, rows = vars(id_model))
p <- p + facet_grid(id_model~Treatment)
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p
## write output
file2write <- paste0(dir_out, "umap.", "allclusters.", "bySample.", id_integration, ".png")
png(filename = file2write, width = 2200, height = 1000, res = 150)
print(p)
dev.off()

# make umap by sample with all clusters together ---------------------------------------------------------------
## make plot data
plotdata_df <- barcode2cluster_df %>%
  mutate(Id_Manual_Cluster = factor(Id_Manual_Cluster)) %>%
  mutate(id_model = str_split_fixed(string = Id_Sample, pattern = "-", n = 3)[,1]) %>%
  mutate(treatment = str_split_fixed(string = Id_Sample, pattern = "-", n = 3)[,3])
unique(plotdata_df$treatment)
## make plot
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = Id_Manual_Cluster), 
                    size = 0.2, alpha = 0.8)
p <- p + scale_color_manual(values = colors_cluster)
p <- p + facet_grid(rows = vars(id_model), cols = vars(treatment))
# p <- p + facet_grid(cols = vars(id_model))
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p
## write output
file2write <- paste0(dir_out, "umap.", "allclusters.", "bySample.", id_integration, ".png")
png(filename = file2write, width = 4200, height = 2000, res = 150)
print(p)
dev.off()


