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
wb_value_df = read_xlsx("./Resources/Western_Blot/Image J quant_05312023_transposed.xlsx")

# set parameter -----------------------------------------------------------
phospho_plot = "phospho_4EBP1"
total_plot = "total_4EBP1"
loadcontrol_plot = "loadcontrol"
stderror <- function(x) sd(x)/sqrt(length(x))

# plot --------------------------------------------------------------------
wb_value_filtered_df = wb_value_df[!is.na(wb_value_df[,phospho_plot]) & !is.na(wb_value_df[,total_plot]),]
wb_value_filtered_df[,"phospho"] = wb_value_filtered_df[,phospho_plot]
wb_value_filtered_df[,"total_protein"] = wb_value_filtered_df[,total_plot]
wb_value_filtered_df = wb_value_filtered_df %>%
  mutate(phospho.bytotal = phospho/total_protein) %>%
  mutate(phospho.bytotal.bylc = (phospho/total_protein)/loadcontrol)

phospho.bymodel_vec = NULL
phospho.bytotal.bymodel_vec = NULL
phospho.bytotal.bylc.bymodel_vec = NULL
for (date_tmp in unique(wb_value_filtered_df$Date)) {
  wb_value_tmp_df = wb_value_filtered_df[wb_value_filtered_df$Date == date_tmp,]
  phospho.bymodel_vec = c(phospho.bymodel_vec, 100*(wb_value_tmp_df$phospho)/wb_value_tmp_df$phospho[(wb_value_tmp_df$Model == "RESL5") & (wb_value_tmp_df$Treatment_group == "Control")])
  phospho.bytotal.bymodel_vec = c(phospho.bytotal.bymodel_vec, 100*(wb_value_tmp_df$phospho.bytotal)/wb_value_tmp_df$phospho.bytotal[(wb_value_tmp_df$Model == "RESL5") & (wb_value_tmp_df$Treatment_group == "Control")])
  phospho.bytotal.bylc.bymodel_vec = c(phospho.bytotal.bylc.bymodel_vec, 100*(wb_value_tmp_df$phospho.bytotal.bylc)/wb_value_tmp_df$phospho.bytotal.bylc[(wb_value_tmp_df$Model == "RESL5") & (wb_value_tmp_df$Treatment_group == "Control")])
}
wb_value_filtered_df$phospho.bymodel = phospho.bymodel_vec
wb_value_filtered_df$phospho.bytotal.bymodel = phospho.bytotal.bymodel_vec
wb_value_filtered_df$phospho.bytotal.bylc.bymodel = phospho.bytotal.bylc.bymodel_vec

plot_data_df = wb_value_filtered_df %>%
  group_by(Model, Treatment_group) %>%
  summarise(phospho.mean = mean(phospho.bymodel), 
            phospho.sem = stderror(phospho.bymodel), 
            phospho.bytotal.mean = mean(phospho.bytotal.bymodel),
            phospho.bytotal.sem = stderror(phospho.bytotal.bymodel), 
            phospho.bytotal.bylc.mean = mean(phospho.bytotal.bylc.bymodel),
            phospho.bytotal.bylc.sem = stderror(phospho.bytotal.bylc.bymodel)) %>%
  mutate(x_plot = paste0(Model, "_", Treatment_group))
plot_data_df$Treatment_group = factor(plot_data_df$Treatment_group, levels = c("Control", "Cabozantinib", "Sapanisertib", "Combo"))

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
pdf(file2write, width = 5, height = 3.5, useDingbats = F)
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
p <- p + ylab(paste0(phospho_plot, "\nnormalized by total protein", "\nnormalized by RESL5 control"))
p <- p + theme(axis.text.x = element_blank(), 
               axis.title.x = element_blank(),
               axis.ticks.x = element_blank())
pdf.options(reset = TRUE, onefile = FALSE)
file2write <- paste0(dir_out, phospho_plot, ".phospho.bytotal",  ".pdf")
pdf(file2write, width = 5, height = 3.5, useDingbats = F)
print(p)
dev.off()

# plot phospho normalized by total normalized by loading control ------------------------------------------------------------
p <- ggplot(data = plot_data_df, aes(x = Treatment_group, y = phospho.bytotal.bylc.mean))
p <- p + geom_bar(aes(fill = Treatment_group), 
                  stat = "identity")
p <- p + geom_errorbar(aes(ymin=phospho.bytotal.bylc.mean-phospho.bytotal.bylc.sem, 
                           ymax=phospho.bytotal.bylc.mean+phospho.bytotal.bylc.sem), 
                       width=0.5, alpha = 0.5, size = 0.5)
p <- p + facet_grid(.~Model, scales = "free_x")
p <- p + ylab(paste0(phospho_plot, "\nnormalized by total protein, loading control, ", "\n and RESL5 control"))
p <- p + theme(axis.text.x = element_blank(), 
               axis.title.x = element_blank(),
               axis.ticks.x = element_blank())
pdf.options(reset = TRUE, onefile = FALSE)
file2write <- paste0(dir_out, phospho_plot, ".phospho.bytotal.byloadingcontrol",  ".pdf")
pdf(file2write, width = 5, height = 3.5, useDingbats = F)
print(p)
dev.off()

