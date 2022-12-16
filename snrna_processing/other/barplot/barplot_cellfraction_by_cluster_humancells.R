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
  mutate(frac_bysample_bycluster = 100*(count_bysample_bycluster/count_bysample)) %>%
  mutate(id_model = str_split_fixed(string = orig.ident, pattern = "-", n = 3)[,1]) %>%
  mutate(id_model = ifelse(id_model == "RESL10F", "RESL10", "RESL5")) %>%
  mutate(treatment_group = str_split_fixed(string = orig.ident, pattern = "-", n = 3)[,3]) %>%
  mutate(treatment_group = gsub(pattern = "2", replacement = "", x = treatment_group)) %>%
  mutate(Treatment = ifelse(treatment_group == "CT", "Control",
                            ifelse(treatment_group == "Cabo", "Cabozantinib",
                                   ifelse(treatment_group == "Sap", "Sapanisertib", "Cabozantinib+\nSapanisertib")))) %>%
  mutate(Treatment = factor(x = Treatment, levels = c("Control", "Cabozantinib+\nSapanisertib", "Cabozantinib", "Sapanisertib"))) %>%
  mutate(cluster = paste0("MC", (seurat_clusters+1))) %>%
  mutate(cluster = factor(cluster))

## colors for treatment
color_red <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[1]
color_green <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[3]
color_yellow <- RColorBrewer::brewer.pal(n = 7, name = "Set2")[6]
color_grey <- "grey50"
colors_treatment <- c("Control" = color_grey, "Cabozantinib" = color_red, "Sapanisertib" = color_green, "Cabozantinib+\nSapanisertib" = color_yellow)
## define size
fontsize_plot <- 12

# make plot for CT-Sap-Combo ---------------------------------------------------------
for (cluster_tmp in unique(plot_data_df$cluster)) {
  p <- ggplot(data = subset(plot_data_df, cluster == cluster_tmp), mapping = aes(x = Treatment, y = frac_bysample_bycluster, fill = Treatment))
  p <- p + geom_bar(stat = "identity")
  # p <- p + geom_point()
  p <- p + scale_fill_manual(values = colors_treatment)
  p <- p + facet_grid(cols = vars(id_model), scales = "free")
  p <- p + ylab(label = paste0("% of ", cluster_tmp, " in the tumor cells"))
  p <- p + theme_classic()
  p <- p + guides(colour = guide_legend(ncol = 1))
  p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = fontsize_plot),
                 axis.text.x = element_blank(),
                 strip.background = element_rect(color = NA), strip.text = element_text(size = fontsize_plot),
                 axis.text.y = element_text(color = "black", size = fontsize_plot), 
                 legend.text = element_text(size = fontsize_plot), legend.title = element_text(size = fontsize_plot),
                 axis.ticks.x = element_blank(), 
                 axis.line.x = element_blank())
  # file2write <- paste0(dir_out, cluster_tmp, ".png")
  # png(filename = file2write, width = 1000, height = 600, res = 150)
  # print(p)
  # dev.off()
  file2write <- paste0(dir_out, cluster_tmp, ".pdf")
  pdf(file2write, width = 4, height = 2.5, useDingbats = F)
  print(p)
  dev.off()
}
