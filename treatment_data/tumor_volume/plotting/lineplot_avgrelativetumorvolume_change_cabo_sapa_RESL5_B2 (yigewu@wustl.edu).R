# Yige Wu @WashU Jun 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
source("./ccRCC_drug_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies --------------------------------------------------
## input the treatment groups
treatmentgroups_process <- c("CT", "Cab", "Sap", "Cab+Sap")
## set sample id
id_sample <- "RESL5_B2"
## input formatted measurement data
avg_relative_volume_df <- fread(data.table = F, input = "./Resources/Analysis_Results/treatment_data/tumor_volume/process_tumor_volume/20210128.v1/RESL5_B2_XY_20210128/RESL5_B2_XY_20210128.RelativeTumoreVolume.Average.ByTreatmentGroup.20210128.v1.tsv")

# make plot data ----------------------------------------------------------
avg_relative_volume_df <- avg_relative_volume_df %>%
  mutate(Date = as.Date(Date, "%Y-%m-%d")) %>%
  filter(Treatment_Status == "Treatment on")
avg_relative_volume_df$Treatment_Days = (avg_relative_volume_df$Date - avg_relative_volume_df$Date[1])
plot_data_wide_df <- avg_relative_volume_df[, c("Treatment_Days", treatmentgroups_process)]
plot_data_long_df <- reshape2::melt(data = plot_data_wide_df, id.var = c("Treatment_Days"))
plot_data_long_df <- plot_data_long_df %>%
  rename(Avg_Relative_Volume = value) %>%
  rename(Treatment_Group = variable) %>%
  mutate(Avg_Relative_Volume_Change = (Avg_Relative_Volume - 100))
  

# specify plotting parameters ---------------------------------------------
color_red <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[1]
color_green <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[3]
color_yellow <- RColorBrewer::brewer.pal(n = 7, name = "Set2")[6]
color_grey <- "grey50"

# plot complete plot --------------------------------------------------------------------
p <- ggplot(data = plot_data_long_df, mapping = aes(x = Treatment_Days, y = Avg_Relative_Volume_Change, group = Treatment_Group, color = Treatment_Group))
p <- p + geom_line()
p <- p + geom_point()
p <- p + scale_color_manual(values = c("CT" = color_grey, "Cab" = color_red, "Sap" = color_green, "Cab+Sap" = color_yellow))
# p <- p + coord_cartesian(ylim = c(-41, 100), xlim = c(0, 60), expand = FALSE)
# p <- p + scale_y_continuous(breaks = seq(-40, 100, 20))
# p <- p + scale_x_continuous(expand = c(0,0), breaks = seq(0, 60, 10))
p <- p + theme_classic(base_size = 16)
p <- p + geom_hline(aes(yintercept = 0), linetype = 2)
p <- p + theme(panel.grid.major.y = element_line(colour = "grey90"))
p <- p + ylab("Tumor volume change (%)") + xlab("Treatment time (days)")
p
file2write <- paste0(dir_out, id_sample, ".", "complete.", "png")
png(filename = file2write, width = 1000, height = 1000, res = 150)
print(p)
dev.off()

# plot 0 - 100 --------------------------------------------------------------------
p <- ggplot(data = plot_data_long_df, mapping = aes(x = Treatment_Days, y = Avg_Relative_Volume_Change, group = Treatment_Group, color = Treatment_Group))
p <- p + geom_line()
p <- p + geom_point()
p <- p + scale_color_manual(values = c("CT" = color_grey, "Cab" = color_red, "Sap" = color_green, "Cab+Sap" = color_yellow))
p <- p + coord_cartesian(ylim = c(-41, 100), xlim = c(0, 60), expand = FALSE)
p <- p + scale_y_continuous(breaks = seq(-40, 100, 20))
p <- p + scale_x_continuous(expand = c(0,0), breaks = seq(0, 60, 10))
p <- p + theme_classic(base_size = 16)
p <- p + geom_hline(aes(yintercept = 0), linetype = 2)
p <- p + theme(panel.grid.major.y = element_line(colour = "grey90"))
p <- p + ylab("Tumor volume change (%)") + xlab("Treatment time (days)")
p
file2write <- paste0(dir_out, id_sample, ".", "y_below100.", "png")
png(filename = file2write, width = 800, height = 600, res = 150)
print(p)
dev.off()

# plot > 100 --------------------------------------------------------------------
p <- ggplot(data = plot_data_long_df, mapping = aes(x = Treatment_Days, y = Avg_Relative_Volume_Change, group = Treatment_Group, color = Treatment_Group))
p <- p + geom_line()
p <- p + geom_point()
p <- p + scale_color_manual(values = c("CT" = color_grey, "Cab" = color_red, "Sap" = color_green, "Cab+Sap" = color_yellow))
p <- p + coord_cartesian(ylim = c(100, 700), xlim = c(0, 60), expand = FALSE)
p <- p + scale_y_continuous(breaks = seq(100, 700, 200))
p <- p + scale_x_continuous(expand = c(0,0), breaks = seq(0, 60, 10))
p <- p + theme_classic(base_size = 16)
p <- p + geom_hline(aes(yintercept = 0), linetype = 2)
p <- p + theme(panel.grid.major.y = element_line(colour = "grey90"),
               axis.text.x = element_blank(), axis.ticks.x = element_blank(),
               axis.line.x = element_blank())
p <- p + ylab("Tumor volume change (%)")
p <- p + theme(axis.title.x = element_blank())
p
file2write <- paste0(dir_out, id_sample, ".", "y_over100.", "png")
png(filename = file2write, width = 800, height = 200, res = 150)
print(p)
dev.off()

