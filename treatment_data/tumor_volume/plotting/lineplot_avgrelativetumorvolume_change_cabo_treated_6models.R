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
avg_relative_volume_df <- read_excel(path = "./Resources/Treatment_Lists/Treatment_Curves/CaboSap_TreatmentCurve_Data.121321.xlsx")

# make plot data ----------------------------------------------------------
plot_data_long_df <- avg_relative_volume_df %>%
  filter(Treatment == "Cabozantinib") %>%
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
colors_plot <- RColorBrewer::brewer.pal(n = 6, name = "Dark2")
colors_plot <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[-6]
names(colors_plot) <- c("RESL5", "RESL10", "RESL12", "RESL4", "RESL11", "RESL3")

# plot complete plot --------------------------------------------------------------------
p <- ggplot(data = plot_data_long_df, mapping = aes(x = Treatment_Days, y = Avg_Relative_Volume_Change, group = Model, color = Model))
p <- p + geom_line()
p <- p + geom_point()
p <- p + geom_errorbar(aes(ymin=Avg_Relative_Volume_Change-STDEV_relative_tumor_volumn, 
                           ymax=Avg_Relative_Volume_Change+STDEV_relative_tumor_volumn), 
                       width=0.75, alpha = 0.2, size = 0.5)
p <- p + scale_color_manual(values = colors_plot)
p <- p + scale_y_continuous(breaks = seq(-100, 700, 100))
p <- p + scale_x_continuous(expand = c(0,0), breaks = seq(0, 70, 10), limits = c(0, 70))
p <- p + theme_classic(base_size = 16)
p <- p + geom_hline(aes(yintercept = 0), linetype = 2)
p <- p + theme(panel.grid.major.y = element_line(colour = "grey90"))
p <- p + ylab("Tumor volume change (%)") + xlab("Treatment time (days)")
p
file2write <- paste0(dir_out, "Cabozantinib.complete.", "png")
png(filename = file2write, width = 1000, height = 600, res = 150)
print(p)
dev.off()

pdf.options(reset = TRUE, onefile = FALSE)
file2write <- paste0(dir_out, "Cabozantinib.complete.", "pdf")
pdf(file2write, width = 5.5, height = 3, useDingbats = F)
print(p)
dev.off()

# plot complete plot - x axis cut --------------------------------------------------------------------
p <- ggplot(data = plot_data_long_df, mapping = aes(x = Treatment_Days, y = Avg_Relative_Volume_Change, group = Model, color = Model))
p <- p + geom_line()
p <- p + geom_point()
p <- p + geom_errorbar(aes(ymin=Avg_Relative_Volume_Change-STDEV_relative_tumor_volumn, 
                           ymax=Avg_Relative_Volume_Change+STDEV_relative_tumor_volumn), 
                       width=0.75, alpha = 0.2, size = 0.5)
p <- p + scale_color_manual(values = colors_plot)
p <- p + scale_y_continuous(breaks = seq(-100, 700, 100))
p <- p + scale_x_continuous(expand = c(0,0), breaks = seq(0, 70, 10), limits = c(0, 49))
p <- p + theme_classic(base_size = 16)
p <- p + geom_hline(aes(yintercept = 0), linetype = 2)
p <- p + theme(panel.grid.major.y = element_line(colour = "grey90"))
p <- p + ylab("Tumor volume change (%)") + xlab("Treatment time (days)")
p
file2write <- paste0(dir_out, "Cabozantinib.xlimit50.", "png")
png(filename = file2write, width = 1000, height = 600, res = 150)
print(p)
dev.off()

pdf.options(reset = TRUE, onefile = FALSE)
file2write <- paste0(dir_out, "Cabozantinib.xlimit50.", "pdf")
pdf(file2write, width = 5, height = 3, useDingbats = F)
print(p)
dev.off()


# plot complete plot - x axis cut --------------------------------------------------------------------
p <- ggplot(data = plot_data_long_df, mapping = aes(x = Treatment_Days, y = Avg_Relative_Volume_Change, group = Model, color = Model))
p <- p + geom_line()
p <- p + geom_point()
p <- p + geom_errorbar(aes(ymin=Avg_Relative_Volume_Change-STDEV_relative_tumor_volumn, 
                           ymax=Avg_Relative_Volume_Change+STDEV_relative_tumor_volumn), 
                       width=0.75, alpha = 0.2, size = 0.5)
p <- p + scale_color_manual(values = colors_plot)
p <- p + scale_y_continuous(breaks = seq(-100, 700, 100), limits = c(-100, 200))
p <- p + scale_x_continuous(expand = c(0,0), breaks = seq(0, 70, 10), limits = c(0, 32))
p <- p + theme_classic(base_size = 16)
p <- p + geom_hline(aes(yintercept = 0), linetype = 2)
p <- p + theme(panel.grid.major.y = element_line(colour = "grey90"))
p <- p + ylab("Tumor volume change (%)") + xlab("Treatment time (days)")
p
file2write <- paste0(dir_out, "Cabozantinib.xlimit30.", "png")
png(filename = file2write, width = 1000, height = 600, res = 150)
print(p)
dev.off()

pdf.options(reset = TRUE, onefile = FALSE)
file2write <- paste0(dir_out, "Cabozantinib.xlimit30.", "pdf")
pdf(file2write, width = 5, height = 3, useDingbats = F)
print(p)
dev.off()




