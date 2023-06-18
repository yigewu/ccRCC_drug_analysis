# Yige Wu @WashU March 2022

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
packages = c(
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "readxl",
  "ggplot2",
  "ggbreak"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg_name_tmp)
  }
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_drug_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies --------------------------------------------------
## input formatted measurement data
# avg_relative_volume_df <- fread(data.table = F, input = "./Resources/Analysis_Results/treatment_data/tumor_volume/process_tumor_volume/20210128.v1/RESL5_B2_XY_20210128/RESL5_B2_XY_20210128.RelativeTumoreVolume.Average.ByTreatmentGroup.20210128.v1.tsv")
avg_relative_volume_df <- read_excel(path = "./Resources/Treatment_Lists/Treatment_Curves/CaboSap_TreatmentCurve_Data.xlsx")
## set sample id
id_model <- "RESL10"; id_batch <- 1

# make plot data ----------------------------------------------------------
plot_data_long_df <- avg_relative_volume_df %>%
  filter(Model == id_model) %>%
  filter(Batch == id_batch) %>%
  mutate(Treatment_Days = Treatment_days) %>%
  mutate(Avg_Relative_Volume = Relative_tumor_volumn) %>%
  mutate(Treatment_Group = Treatment) %>%
  mutate(Avg_Relative_Volume_Change = (Avg_Relative_Volume - 100))

# specify plotting parameters ---------------------------------------------
color_red <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[1]
color_green <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[3]
color_yellow <- RColorBrewer::brewer.pal(n = 7, name = "Set2")[6]
color_grey <- "grey50"
colors_plot <- c("Control" = color_grey, "Cabozantinib" = color_red, "Sapanisertib" = color_green, "Cabozantinib+Sapanisertib" = color_yellow)

# plot complete plot --------------------------------------------------------------------
p <- ggplot(data = plot_data_long_df, mapping = aes(x = Treatment_Days, y = Avg_Relative_Volume, group = Treatment_Group, color = Treatment_Group))
p <- p + geom_line()
p <- p + geom_point(alpha = 0.7, shape = 16)
p <- p + geom_errorbar(aes(ymin=Avg_Relative_Volume-STDEV_relative_tumor_volumn, 
                           ymax=Avg_Relative_Volume+STDEV_relative_tumor_volumn), 
                       width=1, alpha = 0.7, size = 0.5)
p <- p + scale_color_manual(values = colors_plot)
p <- p + scale_y_continuous(expand = c(0,0), breaks = seq(0, 1000, 100), limits = c(0, 1000))
p <- p + theme_classic(base_size = 16)
p <- p + theme(panel.grid.major.y = element_line(colour = "grey90"), legend.position = "none")
p <- p + scale_x_continuous(expand = c(0,0), breaks = seq(0, 50, 10), limits = c(0, 51))
p <- p + geom_hline(aes(yintercept = 100), linetype = 2)
p <- p + ylab("Relative tumor volume (%)") + xlab("Treatment time (days)")
# p
# file2write <- paste0(dir_out, id_model, ".B", id_batch, ".cab_sap", ".complete.", "png")
# png(filename = file2write, width = 1000, height = 600, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, id_model, ".B", id_batch, ".cab_sap", ".complete.", "pdf")
pdf(file2write, width = 4, height = 3.5, useDingbats = F)
print(p)
dev.off()

# plot complete plot with ggbreak--------------------------------------------------------------------
p <- ggplot(data = plot_data_long_df, mapping = aes(x = Treatment_Days, y = Avg_Relative_Volume, group = Treatment_Group, color = Treatment_Group))
p <- p + geom_line()
p <- p + geom_point()
p <- p + geom_errorbar(aes(ymin=Avg_Relative_Volume-STDEV_relative_tumor_volumn, 
                           ymax=Avg_Relative_Volume+STDEV_relative_tumor_volumn), 
                       width=1, alpha = 0.5, size = 0.5)
p <- p + scale_color_manual(values = colors_plot)
# p <- p + ylim(c(10, 800))
p <- p + scale_x_continuous(expand = c(0,0), breaks = seq(0, 50, 10), limits = c(0, 51))
p <- p + theme_classic(base_size = 16)
# p <- p + scale_y_cut(breaks = 200, which = 1, scales = 1.5)
p <- p + ylim(c(10, 1000))
p <- p + scale_y_cut(breaks = c(200, 600), which = c(1, 2), scales = c(0.5, 1), space = 0.17)

p <- p + theme(panel.grid.major.y = element_line(colour = "grey90"), legend.position = "none")
p <- p + ylab("Relative tumor volume (%)") + xlab("Treatment time (days)")
# p

# file2write <- paste0(dir_out, id_model, ".B", id_batch, ".cab_sap", ".ggbreak.", "png")
# png(filename = file2write, width = 600, height = 600, res = 150)
# print(p)
# dev.off()

pdf.options(reset = TRUE, onefile = FALSE)
file2write <- paste0(dir_out, id_model, ".B", id_batch, ".cab_sap", ".ggbreak.", "pdf")
pdf(file2write, width = 4, height = 3.5, useDingbats = F)
print(p)
dev.off()

# plot complete plot with ggbreak at 100--------------------------------------------------------------------
p <- ggplot(data = plot_data_long_df, mapping = aes(x = Treatment_Days, y = Avg_Relative_Volume, group = Treatment_Group, color = Treatment_Group))
p <- p + geom_line()
p <- p + geom_point(alpha = 0.7, shape = 16)
p <- p + geom_errorbar(aes(ymin=Avg_Relative_Volume-STDEV_relative_tumor_volumn, 
                           ymax=Avg_Relative_Volume+STDEV_relative_tumor_volumn), 
                       width=1, alpha = 0.5, size = 0.5)
p <- p + scale_color_manual(values = colors_plot)
p <- p + ylim(c(10, 1000))
p <- p + scale_x_continuous(expand = c(0,0), breaks = seq(0, 50, 10), limits = c(0, 51))
p <- p + theme_classic(base_size = 16)
p <- p + scale_y_cut(breaks = c(100), scales = 1, space = 0)
p <- p + theme(panel.grid.major.y = element_line(colour = "grey90"), legend.position = "none")
p <- p + ylab("Relative tumor volume (%)") + xlab("Treatment time (days)")
# p
file2write <- paste0(dir_out, id_model, ".B", id_batch, ".cab_sap", ".ggbreak.v2.", "png")
png(filename = file2write, width = 600, height = 600, res = 150)
print(p)
dev.off()

# for (treatment_tmp in c("Cabozantinib", "Sapanisertib", "Cabozantinib+Sapanisertib")) {
#   p <- ggplot(data = subset(plot_data_long_df, Treatment_Group %in% c("Control", treatment_tmp)), mapping = aes(x = Treatment_Days, y = Avg_Relative_Volume, group = Treatment_Group, color = Treatment_Group))
#   p <- p + geom_line()
#   p <- p + geom_point()
#   p <- p + geom_errorbar(aes(ymin=Avg_Relative_Volume-STDEV_relative_tumor_volumn, 
#                              ymax=Avg_Relative_Volume+STDEV_relative_tumor_volumn), 
#                          width=1, alpha = 0.5, size = 0.5)
#   p <- p + scale_color_manual(values = colors_plot)
#   # p <- p + ylim(c(10, 800))
#   p <- p + scale_x_continuous(expand = c(0,0), breaks = seq(0, 50, 10), limits = c(0, 51))
#   p <- p + theme_classic(base_size = 16)
#   # p <- p + scale_y_cut(breaks = 200, which = 1, scales = 1.5)
#   p <- p + ylim(c(10, 1000))
#   p <- p + scale_y_cut(breaks = c(200, 600), which = c(1, 2), scales = c(0.5, 1), space = 0.17)
#   
#   p <- p + theme(panel.grid.major.y = element_line(colour = "grey90"), legend.position = "none")
#   p <- p + ylab("Relative tumor volume (%)") + xlab("Treatment time (days)")
#   # p
#   
#   # file2write <- paste0(dir_out, id_model, ".B", id_batch, ".cab_sap", ".ggbreak.", "png")
#   # png(filename = file2write, width = 600, height = 600, res = 150)
#   # print(p)
#   # dev.off()
#   
#   pdf.options(reset = TRUE, onefile = FALSE)
#   file2write <- paste0(dir_out, id_model, ".B", id_batch, ".", treatment_tmp, ".ggbreak.", "pdf")
#   pdf(file2write, width = 4, height = 3.5, useDingbats = F)
#   print(p)
#   dev.off()
# }
# 
# 
