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
## input median signature scores per cluster
results_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/signature_scores/other/make_median_signaturescores_bygeneset_bycluster_bysample/20220918.v1/median_scores.res05.bycluster.bysample.humancells.8sampleintegrated.20220918.v1.tsv")
View(data.frame(unique(results_df$gene_set)))


# specify gene set to use -------------------------------------------------
geneset_plot <- "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
# geneset_plot <- "REACTOME_CELL_CYCLE"
# geneset_plot <- "KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450"
# geneset_plot <- "KEGG_DRUG_METABOLISM_CYTOCHROME_P450"
geneset_plot <- "HALLMARK_PI3K_AKT_MTOR_SIGNALING"
geneset_plot <- "WP_PI3KAKT_SIGNALING_PATHWAY"
geneset_plot <- "HALLMARK_MTORC1_SIGNALING"
geneset_plot <- "REACTOME_MAPK3_ERK1_ACTIVATION"
geneset_plot <- "REACTOME_MAPK1_ERK2_ACTIVATION"
geneset_plot <- "REACTOME_MET_RECEPTOR_ACTIVATION"
geneset_plot <- "WP_MET_IN_TYPE_1_PAPILLARY_RENAL_CELL_CARCINOMA"

# make plot data ----------------------------------------------------------
plot_data_df <- results_df %>%
  filter(gene_set %in% geneset_plot) %>%
  mutate(cluster = str_split_fixed(string = group, pattern = "\\-", n = 4)[,4]) %>%
  mutate(model = str_split_fixed(string = group, pattern = "\\-", n = 4)[,1]) %>%
  mutate(treatment = str_split_fixed(string = group, pattern = "\\-", n = 4)[,3]) %>%
  mutate(treatment = gsub(x = treatment, pattern = '[0-9]', replacement = "")) %>%
  mutate(treatment_group = treatment)
plot_data_df$model[plot_data_df$model == "RESL5E"] <- "RESL5"
plot_data_df$model[plot_data_df$model == "RESL10F"] <- "RESL10"
plot_data_df$treatment_group[plot_data_df$treatment == "CT"] <- "Control"
plot_data_df$treatment_group[plot_data_df$treatment == "Cabo"] <- "Cabozantinib"
plot_data_df$treatment_group[plot_data_df$treatment == "Sap"] <- "Sapanisertib"
plot_data_df$treatment_group[plot_data_df$treatment == "Cabo_Sap"] <- "Cabozantinib+\nSapanisertib"
plot_data_df$treatment_group <- factor(x = plot_data_df$treatment_group, levels = c("Control", "Cabozantinib+\nSapanisertib", "Cabozantinib", "Sapanisertib"))
plot_data_df <- plot_data_df %>%
  mutate(cluster = as.numeric(cluster)) %>%
  mutate(cluster_name = paste0("MC", (cluster+1)))

## colors for treatment
color_red <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[1]
color_green <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[3]
color_yellow <- RColorBrewer::brewer.pal(n = 7, name = "Set2")[6]
color_grey <- "grey50"
colors_treatment <- c("Control" = color_grey, "Cabozantinib" = color_red, "Sapanisertib" = color_green, "Cabozantinib+\nSapanisertib" = color_yellow)
## colors for clusters
clusters_uniq_vec <- unique(plot_data_df$cluster_name)
colors_bycluster <- Polychrome::palette36.colors(n = (length(clusters_uniq_vec)+1))[-2]
names(colors_bycluster) <- clusters_uniq_vec
## define size
fontsize_plot <- 12

# make boxplot for all treatment groups ---------------------------------------------------------
p <- ggplot()
p <- p + geom_boxplot(data = plot_data_df, mapping = aes(x = treatment_group, y = value, fill = treatment_group), alpha = 0.5)
p <- p + geom_point(data = plot_data_df, mapping = aes(x = treatment_group, y = value), color = "black", shape = 16)
p <- p + ylab(paste0(geneset_plot, " score"))
p <- p + scale_fill_manual(values = colors_treatment)
p <- p + facet_grid(cols = vars(model), scales = "free")
p <- p + theme_classic()
# p <- p + guides(colour = guide_legend(ncol = 1))
p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = fontsize_plot),
               axis.text.x = element_blank(),
               strip.background = element_rect(color = NA), strip.text = element_text(size = fontsize_plot),
               axis.text.y = element_text(color = "black", size = fontsize_plot),
               legend.text = element_text(size = fontsize_plot), legend.title = element_text(size = fontsize_plot),
               axis.ticks.x = element_blank())
p
file2write <- paste0(dir_out, geneset_plot, ".boxplot.pdf")
pdf(file2write, width = 5, height = 2.5, useDingbats = F)
print(p)
dev.off()

# make boxplot for CT + combo ---------------------------------------------------------
library(ggpubr)

p <- ggplot(data = subset(plot_data_df, treatment_group %in% c("Control", "Cabozantinib+\nSapanisertib")))
p <- p + geom_boxplot( mapping = aes(x = treatment_group, y = value, fill = treatment_group))
# p <- p + geom_point( mapping = aes(x = treatment_group, y = value), color = "black")
p <- p + geom_dotplot(mapping = aes(x = treatment_group, y = value), color = NA, fill = "black", binaxis='y', stackdir='center', dotsize=1, alpha = 0.5, shape = 16)
p <- p + stat_compare_means(mapping = aes(x = treatment_group, y = value, label = paste0("p =", ..p.format..)), method = "wilcox")
# p <- p + ylab(paste0(geneset_plot, " score"))
p <- p + ylab(paste0("EMT score by cluster"))
p <- p + scale_fill_manual(values = colors_treatment[c("Control", "Cabozantinib+\nSapanisertib")])
p <- p + facet_grid(cols = vars(model), scales = "free")
p <- p + ylim(c(0.35, 1.05))
p <- p + theme_classic()
# p <- p + guides(colour = guide_legend(ncol = 1))
p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = fontsize_plot),
               axis.text.x = element_blank(),
               strip.background = element_rect(color = NA), strip.text = element_text(size = fontsize_plot),
               axis.text.y = element_text(color = "black", size = fontsize_plot),
               legend.text = element_text(size = fontsize_plot), legend.title = element_text(size = fontsize_plot),
               axis.ticks.x = element_blank())
file2write <- paste0(dir_out, geneset_plot, ".boxplot.selected.pdf")
pdf(file2write, width = 4, height = 2.5, useDingbats = F)
print(p)
dev.off()

# make boxplot for all treatment groups ---------------------------------------------------------
p <- ggplot(data = subset(plot_data_df, treatment_group %in% c("Control", "Cabozantinib+\nSapanisertib")))
p <- p + geom_point(mapping = aes(x = treatment_group, y = value, color = cluster_name))
p <- p + geom_line(mapping = aes(x = treatment_group, y = value, color = cluster_name, group = cluster_name))
p <- p + ylab(paste0(geneset_plot, " score"))
p <- p + scale_color_manual(values = colors_bycluster)
p <- p + facet_grid(cols = vars(model), scales = "free")
p <- p + theme_classic()
# p <- p + guides(colour = guide_legend(ncol = 1))
p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = fontsize_plot),
               axis.text.x = element_blank(),
               strip.background = element_rect(color = NA), strip.text = element_text(size = fontsize_plot),
               axis.text.y = element_text(color = "black", size = fontsize_plot),
               legend.text = element_text(size = fontsize_plot), legend.title = element_text(size = fontsize_plot),
               axis.ticks.x = element_blank())
p
file2write <- paste0(dir_out, geneset_plot, ".lineplot.pdf")
pdf(file2write, width = 5, height = 2.5, useDingbats = F)
print(p)
dev.off()
