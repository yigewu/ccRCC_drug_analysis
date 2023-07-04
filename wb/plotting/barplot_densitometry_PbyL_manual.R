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
  "ggpubr",
  "rstatix"
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
wb_value_df = read_xlsx("./Resources/Western_Blot/Image J quant_05312023_transposed.xlsx")

# set parameter -----------------------------------------------------------
stderror <- function(x) sd(x)/sqrt(length(x))
color_red <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[1]
color_green <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[3]
color_yellow <- RColorBrewer::brewer.pal(n = 7, name = "Set2")[6]
color_grey <- "grey50"
colors_plot <- c("Control" = color_grey, "Cabozantinib" = color_red, "Sapanisertib" = color_green, "Combo" = color_yellow)

# plot --------------------------------------------------------------------
# phospho_total_plot_df = data.frame(phospho = c("phospho_AKT_473", "phospho_AKT_308", "phospho_ERK", "phospho_4EBP1", "phospho_RPS6"),
#                                    total = c("total_AKT", "total_AKT", "total_ERK", "total_4EBP1", "total_RPS6"))
wb_value_df = wb_value_df %>%
  filter(Date != "060723")
phospho_total_plot_df = data.frame(phospho = c("phospho_AKT_473", "phospho_AKT_308", "phospho_ERK", "phospho_RPS6"),
                                   total = c("total_AKT", "total_AKT", "total_ERK", "total_RPS6"))

# do phospho-AKT-473 ------------------------------------------------------
i = 1
phospho_plot = phospho_total_plot_df$phospho[i]
total_plot = phospho_total_plot_df$total[i]
wb_value_filtered_df = wb_value_df[!is.na(wb_value_df[,phospho_plot]) & !is.na(wb_value_df[,total_plot]),]
wb_value_filtered_df[,"phospho"] = wb_value_filtered_df[,phospho_plot]
wb_value_filtered_df[,"total_protein"] = wb_value_filtered_df[,total_plot]
wb_value_filtered_df = wb_value_filtered_df %>%
  mutate(phospho.bytotal = (phospho/total_protein)*100) %>%
  mutate(phospho.bylc = (phospho/loadcontrol)*100) %>%
  mutate(phospho.bytotal.bylc = (phospho.bytotal/loadcontrol)*100)

plot_data_df = wb_value_filtered_df
plot_data_df$Treatment_group = factor(x = plot_data_df$Treatment_group, levels = c("Control", "Cabozantinib", "Sapanisertib", "Combo"))
plot_data_df$Model = factor(plot_data_df$Model, levels = c("RESL4", "RESL10", "RESL5", "RESL12", "RESL3", "RESL11"))

stat.test <- plot_data_df %>%
  group_by(Model) %>%
  t_test(phospho.bylc ~ Treatment_group, paired = T, 
         comparisons = list(c("Control", "Cabozantinib"), 
                            c("Control", "Sapanisertib"),
                            c("Control", "Combo"),
                            c("Cabozantinib", "Combo"),
                            c("Sapanisertib", "Combo"))) %>%
  adjust_pvalue() %>%
  add_significance("p") %>%
  filter(p < 0.1)
stat.test <- stat.test %>% add_xy_position(fun = "mean_se", x = "Treatment_group")
stat.test$y.position = c(180, 180, 100, 120, 140, 30, 50, 70, 30)

p <- ggbarplot(data = plot_data_df,
               x = "Treatment_group", y = "phospho.bylc", fill = "Treatment_group",
               facet.by = "Model", nrow = 1,
               add = "mean_se", position = position_dodge())
p <- p + ylim(c(0, 220))
p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01)
p <- p + scale_fill_manual(values = colors_plot)
p <- p + ylab(paste0(phospho_plot, " normalized by\nloading control"))
p <- p + theme(axis.text.x = element_blank(),
               axis.title.x = element_blank(),
               axis.ticks.x = element_blank())
pdf.options(reset = TRUE, onefile = FALSE)
file2write <- paste0(dir_out, phospho_plot, ".phospho.bylc.ttest",  ".pdf")
pdf(file2write, width = 7, height = 3.5, useDingbats = F)
print(p)
dev.off()

# do phospho-AKT-308 ------------------------------------------------------
i = 2
phospho_plot = phospho_total_plot_df$phospho[i]
total_plot = phospho_total_plot_df$total[i]
wb_value_filtered_df = wb_value_df[!is.na(wb_value_df[,phospho_plot]) & !is.na(wb_value_df[,total_plot]),]
wb_value_filtered_df[,"phospho"] = wb_value_filtered_df[,phospho_plot]
wb_value_filtered_df[,"total_protein"] = wb_value_filtered_df[,total_plot]
wb_value_filtered_df = wb_value_filtered_df %>%
  mutate(phospho.bytotal = (phospho/total_protein)*100) %>%
  mutate(phospho.bylc = (phospho/loadcontrol)*100) %>%
  mutate(phospho.bytotal.bylc = (phospho.bytotal/loadcontrol)*100)

plot_data_df = wb_value_filtered_df
plot_data_df$Treatment_group = factor(x = plot_data_df$Treatment_group, levels = c("Control", "Cabozantinib", "Sapanisertib", "Combo"))
plot_data_df$Model = factor(plot_data_df$Model, levels = c("RESL4", "RESL10", "RESL5", "RESL12", "RESL3", "RESL11"))

stat.test <- plot_data_df %>%
  group_by(Model) %>%
  t_test(phospho.bylc ~ Treatment_group, paired = T, 
         comparisons = list(c("Control", "Cabozantinib"), 
                            c("Control", "Sapanisertib"),
                            c("Control", "Combo"),
                            c("Cabozantinib", "Combo"),
                            c("Sapanisertib", "Combo"))) %>%
  adjust_pvalue() %>%
  add_significance("p") %>%
  filter(p <= 0.1)
stat.test <- stat.test %>% add_xy_position(fun = "mean_se", x = "Treatment_group")
# stat.test$y.position = c(50, 70, 90, 110, 50, 70, 90)

p <- ggbarplot(data = plot_data_df,
               x = "Treatment_group", y = "phospho.bylc", fill = "Treatment_group",
               facet.by = "Model", nrow = 1,
               add = "mean_se", position = position_dodge())
# p <- p + ylim(c(0, 170))
p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01)
p <- p + scale_fill_manual(values = colors_plot)
p <- p + ylab(paste0(phospho_plot, " normalized by\nloading control"))
# p <- p + ylim(c(0, 6.5))
p <- p + theme(axis.text.x = element_blank(),
               axis.title.x = element_blank(),
               axis.ticks.x = element_blank())
pdf.options(reset = TRUE, onefile = FALSE)
file2write <- paste0(dir_out, phospho_plot, ".phospho.bylc.ttest",  ".pdf")
pdf(file2write, width = 7, height = 3.5, useDingbats = F)
print(p)
dev.off()

# do phospho-ERK ------------------------------------------------------
i = 3
phospho_plot = phospho_total_plot_df$phospho[i]
total_plot = phospho_total_plot_df$total[i]
wb_value_filtered_df = wb_value_df[!is.na(wb_value_df[,phospho_plot]) & !is.na(wb_value_df[,total_plot]),]
wb_value_filtered_df[,"phospho"] = wb_value_filtered_df[,phospho_plot]
wb_value_filtered_df[,"total_protein"] = wb_value_filtered_df[,total_plot]
wb_value_filtered_df = wb_value_filtered_df %>%
  mutate(phospho.bytotal = (phospho/total_protein)*100) %>%
  mutate(phospho.bylc = (phospho/loadcontrol)*100) %>%
  mutate(phospho.bytotal.bylc = (phospho.bytotal/loadcontrol)*100)

plot_data_df = wb_value_filtered_df
plot_data_df$Treatment_group = factor(x = plot_data_df$Treatment_group, levels = c("Control", "Cabozantinib", "Sapanisertib", "Combo"))
plot_data_df$Model = factor(plot_data_df$Model, levels = c("RESL4", "RESL10", "RESL5", "RESL12", "RESL3", "RESL11"))

stat.test <- plot_data_df %>%
  group_by(Model) %>%
  t_test(phospho.bylc ~ Treatment_group, paired = T, 
         comparisons = list(c("Control", "Cabozantinib"), 
                            c("Control", "Sapanisertib"),
                            c("Control", "Combo"),
                            c("Cabozantinib", "Combo"),
                            c("Sapanisertib", "Combo"))) %>%
  adjust_pvalue() %>%
  add_significance("p") %>%
  filter(p <= 0.1)
stat.test <- stat.test %>% add_xy_position(fun = "mean_se", x = "Treatment_group")
stat.test$y.position[stat.test$Model == "RESL4"] = c(130, 160, 190)
stat.test$y.position[stat.test$Model == "RESL10"] = c(150, 175, 200, 225, 250)
stat.test$y.position[stat.test$Model == "RESL5"] = c(130, 155, 180, 205, 230)
stat.test$y.position[stat.test$Model == "RESL12"] = c(70, 100, 130, 160)
stat.test$y.position[stat.test$Model == "RESL3"] = c(50, 80, 110)
stat.test$y.position[stat.test$Model == "RESL11"] = c(50, 80, 110, 140)

p <- ggbarplot(data = plot_data_df,
               x = "Treatment_group", y = "phospho.bylc", fill = "Treatment_group",
               facet.by = "Model", nrow = 1,
               add = "mean_se", position = position_dodge())
p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, size = 5)
p <- p + scale_fill_manual(values = colors_plot)
p <- p + ylab(paste0("Normalized pERK"))
p <- p + guides(fill = guide_legend(title = ""))
p <- p + theme_classic(base_size = 18)
p <- p + ylim(c(0, 260))
p <- p + theme(axis.text.x = element_blank(),
               axis.title.x = element_blank(),
               axis.ticks.x = element_blank(), legend.position = "bottom")
pdf.options(reset = TRUE, onefile = FALSE)
file2write <- paste0(dir_out, phospho_plot, ".phospho.bylc.ttest",  ".pdf")
pdf(file2write, width = 7, height = 3.5, useDingbats = F)
print(p)
dev.off()

# do phospho-RPS6 ------------------------------------------------------
i = 4
phospho_plot = phospho_total_plot_df$phospho[i]
total_plot = phospho_total_plot_df$total[i]
wb_value_filtered_df = wb_value_df[!is.na(wb_value_df[,phospho_plot]) & !is.na(wb_value_df[,total_plot]),]
wb_value_filtered_df[,"phospho"] = wb_value_filtered_df[,phospho_plot]
wb_value_filtered_df[,"total_protein"] = wb_value_filtered_df[,total_plot]
wb_value_filtered_df = wb_value_filtered_df %>%
  mutate(phospho.bytotal = (phospho/total_protein)*100) %>%
  mutate(phospho.bylc = (phospho/loadcontrol)*100) %>%
  mutate(phospho.bytotal.bylc = (phospho.bytotal/loadcontrol)*100)

plot_data_df = wb_value_filtered_df
plot_data_df$Treatment_group = factor(x = plot_data_df$Treatment_group, levels = c("Control", "Cabozantinib", "Sapanisertib", "Combo"))
plot_data_df$Model = factor(plot_data_df$Model, levels = c("RESL4", "RESL10", "RESL5", "RESL12", "RESL3", "RESL11"))

stat.test <- plot_data_df %>%
  group_by(Model) %>%
  t_test(phospho.bylc ~ Treatment_group, paired = T, 
         comparisons = list(c("Control", "Cabozantinib"), 
                            c("Control", "Sapanisertib"),
                            c("Control", "Combo"),
                            c("Cabozantinib", "Combo"),
                            c("Sapanisertib", "Combo"))) %>%
  adjust_pvalue() %>%
  add_significance("p") %>%
  filter(p <= 0.1)
stat.test <- stat.test %>% add_xy_position(fun = "mean_se", x = "Treatment_group")
stat.test$y.position[stat.test$Model == "RESL4"] = c(130, 150, 170, 190)
stat.test$y.position[stat.test$Model == "RESL10"] = c(130, 150, 170, 190)
stat.test$y.position[stat.test$Model == "RESL5"] = c(130, 150, 170, 190)
stat.test$y.position[stat.test$Model == "RESL12"] = c(130, 150)
stat.test$y.position[stat.test$Model == "RESL3"] = c(100, 120, 140)
stat.test$y.position[stat.test$Model == "RESL11"] = c(100, 120, 140, 160)

p <- ggbarplot(data = plot_data_df,
               x = "Treatment_group", y = "phospho.bylc", fill = "Treatment_group",
               facet.by = "Model", nrow = 1,
               add = "mean_se", position = position_dodge())
p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01)
p <- p + scale_fill_manual(values = colors_plot)
p <- p + ylab(paste0(phospho_plot, " normalized by\nloading control"))
p <- p + ylim(c(0, 200))
p <- p + theme(axis.text.x = element_blank(),
               axis.title.x = element_blank(),
               axis.ticks.x = element_blank())
pdf.options(reset = TRUE, onefile = FALSE)
file2write <- paste0(dir_out, phospho_plot, ".phospho.bylc.ttest",  ".pdf")
pdf(file2write, width = 7, height = 3.5, useDingbats = F)
print(p)
dev.off()
