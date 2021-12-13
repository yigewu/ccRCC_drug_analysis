# Yige Wu @WashU Jun 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
source("./ccRCC_drug_analysis/plotting.R")
library(readxl)
library(ggbreak)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies --------------------------------------------------
## input formatted measurement data
# avg_relative_volume_df <- fread(data.table = F, input = "./Resources/Analysis_Results/treatment_data/tumor_volume/process_tumor_volume/20210128.v1/RESL5_B2_XY_20210128/RESL5_B2_XY_20210128.RelativeTumoreVolume.Average.ByTreatmentGroup.20210128.v1.tsv")
avg_relative_volume_df <- read_excel(path = "./Resources/Treatment_Lists/Treatment_Curves/CaboSap_TreatmentCurve_Data.xlsx")
## set sample id
id_model <- "RESL4"; id_batch <- 1

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
p <- ggplot(data = plot_data_long_df, mapping = aes(x = Treatment_Days, y = Avg_Relative_Volume_Change, group = Treatment_Group, color = Treatment_Group))
p <- p + geom_line()
p <- p + geom_point()
p <- p + geom_errorbar(aes(ymin=Avg_Relative_Volume_Change-STDEV_relative_tumor_volumn, 
                           ymax=Avg_Relative_Volume_Change+STDEV_relative_tumor_volumn), 
                       width=4,
                position=position_dodge(0.05))
p <- p + scale_color_manual(values = colors_plot)
p <- p + scale_y_continuous(breaks = seq(-100, 700, 100))
# p <- p + scale_x_continuous(expand = c(0,0), breaks = seq(0, 70, 10), limits = c(0, 70))
p <- p + theme_classic(base_size = 16)
p <- p + geom_hline(aes(yintercept = 0), linetype = 2)
p <- p + theme(panel.grid.major.y = element_line(colour = "grey90"))
p <- p + ylab("Tumor volume change (%)") + xlab("Treatment time (days)")
# p
file2write <- paste0(dir_out, id_model, ".B", id_batch, ".complete.", "png")
png(filename = file2write, width = 1000, height = 600, res = 150)
print(p)
dev.off()

# plot complete plot with ggbreak--------------------------------------------------------------------
p <- ggplot(data = plot_data_long_df, mapping = aes(x = Treatment_Days, y = Avg_Relative_Volume_Change, group = Treatment_Group, color = Treatment_Group))
p <- p + geom_line()
p <- p + geom_point()
p <- p + geom_errorbar(aes(ymin=Avg_Relative_Volume_Change-STDEV_relative_tumor_volumn, 
                           ymax=Avg_Relative_Volume_Change+STDEV_relative_tumor_volumn), 
                       width=1, alpha = 0.5, size = 0.5)
p <- p + scale_color_manual(values = colors_plot)
# p <- p + scale_y_continuous(breaks = c(seq(-80, 100, 20), seq(100, 700, 200)))
p <- p + ylim(c(-80, 700))
p <- p + scale_x_continuous(expand = c(0,0), breaks = seq(0, 60, 10), limits = c(0, 55))
p <- p + theme_classic(base_size = 16)
p <- p + scale_y_cut(breaks = 100, which = 1, scales = 1.5)
p <- p + theme(panel.grid.major.y = element_line(colour = "grey90"), legend.position = "none")
p <- p + ylab("Tumor volume change (%)") + xlab("Treatment time (days)")
# p

file2write <- paste0(dir_out, id_model, ".B", id_batch, ".ggbreak.", "png")
png(filename = file2write, width = 600, height = 600, res = 150)
print(p)
dev.off()

pdf.options(reset = TRUE, onefile = FALSE)
file2write <- paste0(dir_out, id_model, ".B", id_batch, ".ggbreak.", "pdf")
pdf(file2write, width = 4, height = 3.5, useDingbats = F)
print(p)
dev.off()



