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
  "ggbreak",
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
phospho_plot = "phospho_ERK"
total_plot = "total_ERK"
loadcontrol_plot = "loadcontrol"
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
test_result_df = NULL
for (i in c(1)) {
# for (i in 1:nrow(phospho_total_plot_df)) {
  phospho_plot = phospho_total_plot_df$phospho[i]
  total_plot = phospho_total_plot_df$total[i]
  wb_value_filtered_df = wb_value_df[!is.na(wb_value_df[,phospho_plot]) & !is.na(wb_value_df[,total_plot]),]
  wb_value_filtered_df[,"phospho"] = wb_value_filtered_df[,phospho_plot]
  wb_value_filtered_df[,"total_protein"] = wb_value_filtered_df[,total_plot]
  wb_value_filtered_df = wb_value_filtered_df %>%
    mutate(phospho.bytotal = (phospho/total_protein)*100) %>%
    mutate(phospho.bylc = (phospho/loadcontrol)*100) %>%
    mutate(total.bylc = (total_protein/loadcontrol)*100) %>%
    mutate(phospho.bytotal.bylc = (phospho.bytotal/loadcontrol)*100)
  
  phospho.bymodel_vec = NULL
  phospho.bytotal.bymodel_vec = NULL
  phospho.bytotal.bylc.bymodel_vec = NULL
  for (date_tmp in unique(wb_value_filtered_df$Date)) {
    wb_value_tmp_df = wb_value_filtered_df[wb_value_filtered_df$Date == date_tmp,]
    phospho.bymodel_vec = c(phospho.bymodel_vec, 100*(wb_value_tmp_df$phospho)/wb_value_tmp_df$phospho[(wb_value_tmp_df$Model == "RESL5") & (wb_value_tmp_df$Treatment_group == "Control")])
  }
  wb_value_filtered_df$phospho.bymodel = phospho.bymodel_vec
  wb_value_filtered_df$phospho.bytotal.bylc.bymodel = phospho.bytotal.bylc.bymodel_vec
  
  ## do test
  model_vec = NULL
  treatmentgroup_vec = NULL
  pvalue_vec = NULL
  t_vec = NULL
  fdr_vec = NULL
  for (treatment_group in c("Cabozantinib", "Sapanisertib", "Combo")) {
    pvalue_tmp_vec =  NULL
    for (model_tmp in unique(wb_value_filtered_df$Model)) {
      x = wb_value_filtered_df$phospho.bytotal.bylc[wb_value_filtered_df$Model == model_tmp & wb_value_filtered_df$Treatment_group == treatment_group]
      y = wb_value_filtered_df$phospho.bytotal.bylc[wb_value_filtered_df$Model == model_tmp & wb_value_filtered_df$Treatment_group == "Control"]
      t_result = t.test(x, y, paired = T)
      model_vec = c(model_vec, model_tmp)
      treatmentgroup_vec = c(treatmentgroup_vec, treatment_group)
      pvalue_tmp_vec = c(pvalue_tmp_vec, t_result$p.value)
      t_vec = c(t_vec, t_result$estimate)
    }
    pvalue_vec = c(pvalue_vec, pvalue_tmp_vec)
    fdr_vec = c(fdr_vec, p.adjust(pvalue_tmp_vec, method = "fdr"))
  }
  test_result_tmp_df = data.frame(phospho = phospho_plot, total = total_plot,
                                  model = model_vec, treatment_group = treatmentgroup_vec,
                                  fdr = fdr_vec, p_value = pvalue_vec, t = t_vec)
  test_result_df = rbind(test_result_df, test_result_tmp_df)
  
  plot_data_df = wb_value_filtered_df %>%
    group_by(Model, Treatment_group) %>%
    summarise(phospho.mean = mean(phospho.bymodel), 
              phospho.sem = stderror(phospho.bymodel), 
              phospho.bytotal.mean = mean(phospho.bytotal),
              phospho.bytotal.sem = stderror(phospho.bytotal), 
              total.bylc.mean = mean(total.bylc),
              total.bylc.sem = stderror(total.bylc),
              phospho.bytotal.bylc.mean = mean(phospho.bytotal.bylc),
              phospho.bytotal.bylc.sem = stderror(phospho.bytotal.bylc)) %>%
    mutate(x_plot = paste0(Model, "_", Treatment_group))
  plot_data_df$Treatment_group = factor(plot_data_df$Treatment_group, levels = c("Control", "Cabozantinib", "Sapanisertib", "Combo"))
  plot_data_df$Model = factor(plot_data_df$Model, levels = c("RESL5", "RESL10", "RESL12", "RESL4", "RESL3", "RESL11"))
  
  # plot phospho normalized by total normalized by loading control ------------------------------------------------------------
  # Statistical test
  df = wb_value_filtered_df
  df$Treatment_group = factor(x = df$Treatment_group, levels = c("Control", "Cabozantinib", "Sapanisertib", "Combo"))
  df$Model = factor(df$Model, levels = c("RESL4", "RESL10", "RESL5", "RESL12", "RESL3", "RESL11"))
  
  stat.test <- df %>%
    group_by(Model) %>%
    t_test(phospho.bytotal.bylc ~ Treatment_group, paired = T, 
           comparisons = list(c("Control", "Cabozantinib"), 
                              c("Control", "Sapanisertib"),
                              c("Control", "Combo"),
                              c("Cabozantinib", "Combo"),
                              c("Sapanisertib", "Combo"))) %>%
    adjust_pvalue() %>%
    add_significance("p")
  stat.test <- stat.test %>% add_xy_position(fun = "mean_se", x = "Treatment_group")
  stat.test$p.format <- ifelse(stat.test$p <= 0.1, stat.test$p, "ns")
  
  p <- ggbarplot(data = df,
                 x = "Treatment_group", y = "phospho.bytotal.bylc", fill = "Treatment_group",
                 facet.by = "Model", nrow = 1,
                 add = "mean_se", position = position_dodge())
  p <- p + stat_pvalue_manual(stat.test, label = "p.format", tip.length = 0.01)
  # p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01)
  p <- p + scale_fill_manual(values = colors_plot)
  p <- p + ylab(paste0(phospho_plot, "\nnormalized by total protein, loading control"))
  p <- p + ylim(c(0, max(stat.test$y.position)*1.2))
  p <- p + theme(axis.text.x = element_blank(),
                 axis.title.x = element_blank(),
                 axis.ticks.x = element_blank())
  pdf.options(reset = TRUE, onefile = FALSE)
  file2write <- paste0(dir_out, phospho_plot, ".phospho.bytotal.byloadingcontrol.ttest",  ".pdf")
  pdf(file2write, width = 7, height = 3.5, useDingbats = F)
  print(p)
  dev.off()
  
  stat.test <- df %>%
    group_by(Model) %>%
    wilcox_test(phospho.bytotal.bylc ~ Treatment_group, paired = T, 
           comparisons = list(c("Control", "Cabozantinib"), 
                              c("Control", "Sapanisertib"),
                              c("Control", "Combo"),
                              c("Cabozantinib", "Combo"),
                              c("Sapanisertib", "Combo"))) %>%
    adjust_pvalue() %>%
    add_significance("p")
  stat.test <- stat.test %>% add_xy_position(fun = "mean_se", x = "Treatment_group")
  stat.test$p.format <- ifelse(stat.test$p <= 0.1, stat.test$p, "ns")
  
  p <- ggbarplot(data = df,
                 x = "Treatment_group", y = "phospho.bytotal.bylc", fill = "Treatment_group",
                 facet.by = "Model", nrow = 1,
                 add = "mean_se", position = position_dodge())
  p <- p + stat_pvalue_manual(stat.test, label = "p.format", tip.length = 0.01)
  
  p <- p + scale_fill_manual(values = colors_plot)
  p <- p + ylab(paste0(phospho_plot, "\nnormalized by total protein, loading control"))
  p <- p + theme(axis.text.x = element_blank(),
                 axis.title.x = element_blank(),
                 axis.ticks.x = element_blank())
  pdf.options(reset = TRUE, onefile = FALSE)
  file2write <- paste0(dir_out, phospho_plot, ".phospho.bytotal.byloadingcontrol.wilcox",  ".pdf")
  pdf(file2write, width = 7, height = 3.5, useDingbats = F)
  print(p)
  dev.off()
  
  # plot phospho normalized by total ------------------------------------------------------------
  # Statistical test
  df = wb_value_filtered_df
  df$Treatment_group = factor(x = df$Treatment_group, levels = c("Control", "Cabozantinib", "Sapanisertib", "Combo"))
  df$Model = factor(df$Model, levels = c("RESL4", "RESL10", "RESL5", "RESL12", "RESL3", "RESL11"))
  
  stat.test <- df %>%
    group_by(Model) %>%
    t_test(phospho.bytotal ~ Treatment_group, paired = T, 
           comparisons = list(c("Control", "Cabozantinib"), 
                              c("Control", "Sapanisertib"),
                              c("Control", "Combo"),
                              c("Cabozantinib", "Combo"),
                              c("Sapanisertib", "Combo"))) %>%
    adjust_pvalue() %>%
    add_significance("p")
  stat.test <- stat.test %>% add_xy_position(fun = "mean_se", x = "Treatment_group")
  stat.test$p.format <- ifelse(stat.test$p <= 0.1, stat.test$p, "ns")
  
  p <- ggbarplot(data = df,
                 x = "Treatment_group", y = "phospho.bytotal", fill = "Treatment_group",
                 facet.by = "Model", nrow = 1,
                 add = "mean_se", position = position_dodge())
  p <- p + stat_pvalue_manual(stat.test, label = "p.format", tip.length = 0.01)
  # p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01)
  p <- p + scale_fill_manual(values = colors_plot)
  p <- p + ylab(paste0(phospho_plot, "\nnormalized by total protein"))
  p <- p + ylim(c(0, max(stat.test$y.position)*1.2))
  p <- p + theme(axis.text.x = element_blank(),
                 axis.title.x = element_blank(),
                 axis.ticks.x = element_blank())
  pdf.options(reset = TRUE, onefile = FALSE)
  file2write <- paste0(dir_out, phospho_plot, ".phospho.bytotal.ttest",  ".pdf")
  pdf(file2write, width = 7, height = 3.5, useDingbats = F)
  print(p)
  dev.off()
  
  # plot phospho normalized by loading control ------------------------------------------------------------
  # Statistical test
  df = wb_value_filtered_df
  df$Treatment_group = factor(x = df$Treatment_group, levels = c("Control", "Cabozantinib", "Sapanisertib", "Combo"))
  df$Model = factor(df$Model, levels = c("RESL4", "RESL10", "RESL5", "RESL12", "RESL3", "RESL11"))
  
  stat.test <- df %>%
    group_by(Model) %>%
    t_test(phospho.bylc ~ Treatment_group, paired = T, 
           comparisons = list(c("Control", "Cabozantinib"), 
                              c("Control", "Sapanisertib"),
                              c("Control", "Combo"),
                              c("Cabozantinib", "Combo"),
                              c("Sapanisertib", "Combo"))) %>%
    adjust_pvalue() %>%
    add_significance("p")
  stat.test <- stat.test %>% add_xy_position(fun = "mean_se", x = "Treatment_group")
  stat.test$p.format <- ifelse(stat.test$p <= 0.1, stat.test$p, "ns")
  
  p <- ggbarplot(data = df,
                 x = "Treatment_group", y = "phospho.bylc", fill = "Treatment_group",
                 facet.by = "Model", nrow = 1,
                 add = "mean_se", position = position_dodge())
  p <- p + stat_pvalue_manual(stat.test, label = "p.format", tip.length = 0.01)
  # p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01)
  p <- p + scale_fill_manual(values = colors_plot)
  p <- p + ylab(paste0(phospho_plot, "\nnormalized by loading control"))
  p <- p + ylim(c(0, max(stat.test$y.position)*1.2))
  p <- p + theme(axis.text.x = element_blank(),
                 axis.title.x = element_blank(),
                 axis.ticks.x = element_blank())
  pdf.options(reset = TRUE, onefile = FALSE)
  file2write <- paste0(dir_out, phospho_plot, ".phospho.bylc.ttest",  ".pdf")
  pdf(file2write, width = 7, height = 3.5, useDingbats = F)
  print(p)
  dev.off()
  
  # plot phospho ------------------------------------------------------------
  p <- ggplot(data = plot_data_df, aes(x = Treatment_group, y = phospho.mean))
  p <- p + geom_bar(aes(fill = Treatment_group), 
                    stat = "identity")
  p <- p + geom_errorbar(aes(ymin=phospho.mean-phospho.sem, 
                             ymax=phospho.mean+phospho.sem), 
                         width=0.5, alpha = 0.5, size = 0.5)
  p <- p + facet_grid(.~Model, scales = "free_x")
  p <- p + ylab(paste0(phospho_plot, "\nnormalized by RESL5 control"))
  p <- p + theme(axis.text.x = element_blank(), 
                 axis.title.x = element_blank(),
                 axis.ticks.x = element_blank())
  pdf.options(reset = TRUE, onefile = FALSE)
  file2write <- paste0(dir_out, phospho_plot, ".phospho",  ".pdf")
  pdf(file2write, width = 7, height = 3, useDingbats = F)
  print(p)
  dev.off()
  
  # plot phospho normalized by total ------------------------------------------------------------
  p <- ggplot(data = plot_data_df, aes(x = Treatment_group, y = phospho.bytotal.mean))
  p <- p + geom_bar(aes(fill = Treatment_group), 
                    stat = "identity")
  p <- p + geom_errorbar(aes(ymin=phospho.bytotal.mean-phospho.bytotal.sem, 
                             ymax=phospho.bytotal.mean+phospho.bytotal.sem), 
                         width=0.5, alpha = 0.5, size = 0.5)
  p <- p + facet_grid(.~Model, scales = "free_x")
  p <- p + ylab(paste0(phospho_plot, "\nnormalized by total protein"))
  p <- p + theme(axis.text.x = element_blank(), 
                 axis.title.x = element_blank(),
                 axis.ticks.x = element_blank())
  pdf.options(reset = TRUE, onefile = FALSE)
  file2write <- paste0(dir_out, phospho_plot, ".phospho.bytotal",  ".pdf")
  pdf(file2write, width = 7, height = 3, useDingbats = F)
  print(p)
  dev.off()
  
  # plot phospho normalized by total normalized ------------------------------------------------------------
  # Statistical test
  df = wb_value_filtered_df
  df$Treatment_group = factor(x = df$Treatment_group, levels = c("Control", "Cabozantinib", "Sapanisertib", "Combo"))
  df$Model = factor(df$Model, levels = c("RESL4", "RESL10", "RESL5", "RESL12", "RESL3", "RESL11"))
  
  stat.test <- df %>%
    group_by(Model) %>%
    t_test(phospho.bytotal ~ Treatment_group, paired = T, 
           comparisons = list(c("Control", "Cabozantinib"), 
                              c("Control", "Sapanisertib"),
                              c("Control", "Combo"),
                              c("Cabozantinib", "Combo"),
                              c("Sapanisertib", "Combo"))) %>%
    adjust_pvalue() %>%
    add_significance("p")
  stat.test <- stat.test %>% add_xy_position(fun = "mean_sd", x = "Treatment_group")
  stat.test$p.format <- ifelse(stat.test$p <= 0.1, stat.test$p, "ns")
  
  p <- ggbarplot(data = df,
                 x = "Treatment_group", y = "phospho.bytotal", fill = "Treatment_group",
                 facet.by = "Model", nrow = 1,
                 add = "mean_sd", position = position_dodge())
  p <- p + stat_pvalue_manual(stat.test, label = "p.format", tip.length = 0.01)
  # p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01)
  p <- p + scale_fill_manual(values = colors_plot)
  p <- p + ylab(paste0(phospho_plot, "\nnormalized by total protein, loading control"))
  p <- p + ylim(c(0, max(stat.test$y.position)*1.2))
  p <- p + theme(axis.text.x = element_blank(),
                 axis.title.x = element_blank(),
                 axis.ticks.x = element_blank())
  pdf.options(reset = TRUE, onefile = FALSE)
  file2write <- paste0(dir_out, phospho_plot, ".phospho.bytotal.ttest.mean_sd",  ".pdf")
  pdf(file2write, width = 7, height = 3.5, useDingbats = F)
  print(p)
  dev.off()
  

}
# write.table(x = test_result_df, file = paste0(dir_out, "T_test_result.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)
