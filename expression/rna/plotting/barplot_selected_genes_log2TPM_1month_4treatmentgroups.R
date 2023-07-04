# Yige Wu @ WashU 2023 Jun

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("OmnipathR")
packages = c(
  "plyr",
  "stringr",
  "reshape2",
  "data.table",
  "dplyr"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

# input dependencies ------------------------------------------------------
## input the protein data
exp_df <- fread("./Data_Freeze/v1.dataFreeze.washU_rcc/3.geneExp/v3.20210116/datafreeze.v3.kallisto.geneExp.protein_coding.tsv", data.table = F)
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/0.sample_info/v3.20210116/RCC_PDX_Samples.20210115.v2.xlsx")

# set parameters ----------------------------------------------------------
genes_filter <- c("CCND3", "MYBL2", "FOXM1", "CDK2", "BCL2L12","BCL2L13")
genes_filter <- c("CCND3", "MYBL2", "FOXM1", "BCL2L12")
colnames_id <- colnames(exp_df)[!(colnames(exp_df) %in% sampleinfo_df$Analysis_ID)]

# plot --------------------------------------------------------------------
plot_data_wide_df <- exp_df %>%
  filter(Name %in% genes_filter)
plot_data_long_df <- melt(data = plot_data_wide_df, id.vars = colnames_id)
plot_data_long_df <- as.data.frame(plot_data_long_df)
plot_data_long_df <- merge(x = plot_data_long_df, y = sampleinfo_df, by.x = "variable", by.y = "Analysis_ID", all.x = T)

plot_data_long_df <- plot_data_long_df %>%
  filter(ShortTag %in% c("Control", "Treated.Cabo+Sap", "Treated.Cabo", "Treated.Sap")) %>%
  filter(Treatment.Month == 1) %>%
  mutate(id_model_length = paste0(PairTag, "_", Treatment.Month))
plot_data_long_df$Treatment = mapvalues(plot_data_long_df$ShortTag, 
                                   from = c("Control", "Treated.Cabo", "Treated.Sap", "Treated.Cabo+Sap"),
                                   to = c("Control", "Cabozantinib", "Sapanisertib", "Cabozantinib+Sapanisertib"))
ids_process = unique(plot_data_long_df$id_model_length[plot_data_long_df$Treatment == "Cabozantinib+Sapanisertib"])
plot_data_long_df = plot_data_long_df %>%
  filter(id_model_length %in% ids_process) %>%
  arrange(Name, id_model_length, desc(Treatment))

plot_data_long_df$Treatment = factor(x = plot_data_long_df$Treatment, levels = c("Control", "Cabozantinib", "Sapanisertib", "Cabozantinib+Sapanisertib"))

plot_data_long_df$Model_id <- factor(x =   plot_data_long_df$ModelID, levels = c("RESL5", "RESL10", "RESL11", "RESL3", "RESL12", "RESL4"))
plot_data_long_df$y_plot = log2(plot_data_long_df$value + 1)
# plot_data_long_df$y_plot <- (plot_data_long_df$value - min(plot_data_long_df$value))/(max(plot_data_long_df$value) - min(plot_data_long_df$value))
plot_data_long_df$Name = factor(x = plot_data_long_df$Name, levels = genes_filter)
stat.test <- plot_data_long_df %>%
  group_by(Name) %>%
  t_test(y_plot ~ Treatment, paired = T, 
         comparisons = list(c("Cabozantinib+Sapanisertib", "Control"))) %>%
  adjust_pvalue() %>%
  add_significance("p")
# stat.test <- stat.test %>% add_xy_position(fun = "mean_se", x = "Treatment")
stat.test <- stat.test %>% add_xy_position(fun = "mean_se", x = "Treatment")

## make colors
color_red <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[1]
color_green <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[3]
color_yellow <- RColorBrewer::brewer.pal(n = 7, name = "Set2")[6]
color_grey <- "grey50"
colors_plot <- c("Control" = color_grey, "Cabozantinib" = color_red, "Sapanisertib" = color_green, "Cabozantinib+Sapanisertib" = color_yellow)

# make plot ---------------------------------------------------------------
p <- ggbarplot(data = plot_data_long_df, 
               x = "Treatment", y = "y_plot", fill = "Treatment",
               facet.by = "Name", nrow = 1,
               add = c("mean"),
               position = position_dodge())
p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, label.size = 5)
p <- p + scale_fill_manual(values = colors_plot)
p <- p + labs(y = paste0(" gene expression\n(log2(TPM+1))"))
p <- p + ylim(c(0, 5))
p <- p + theme_classic(base_size = 18)
p <- p + theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.x = element_blank(),
               legend.position = "none")

# write output ------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

## write output
file2write <- paste0(dir_out, "barplot", ".pdf")
pdf(file2write, width = 5, height = 3, useDingbats = F)
print(p)
dev.off()
# file2write <- paste0(dir_out, gene_plot, ".png")
# png(file2write, width = 1500, height = 800, res = 150)
# print(p)
# dev.off()
