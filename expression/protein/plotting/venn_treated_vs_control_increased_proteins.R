# Yige Wu @ WashU 2022 Jan

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
source("./ccRCC_drug_analysis/plotting.R")
## set run id
version_tmp <- 3
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
if (!requireNamespace("eulerr", quietly = TRUE))
  install.packages("eulerr")
library(eulerr)
library(grid)

# input dependencies ------------------------------------------------------
genes_process_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/test/wilcox_paired_diff_protein_treated_vs_control/20220209.v1/Wilcox.Paired.Each_Treated_Group_vs_CT.20220209.v1.tsv")

# make plot data -----------------------------------------------------------
genes_filtered_df <- genes_process_df %>%
  filter(pvalue < 0.05) %>%
  filter(diff_estimate > 0) %>%
  filter(PG.Genes != "") %>%
  unique() %>%
  mutate(treatment_sim = ifelse(group1 == "Cabo", "C",
                                ifelse(group1 == "Sap", "S", "C+S")))

plot_data_df <- dcast(data = genes_filtered_df, formula = PG.Genes ~ treatment_sim, value.var = "diff_estimate")
plot_values_df <- !is.na(plot_data_df[,-1])
plot_data_euler <- euler(plot_values_df)

# plot --------------------------------------------------------------------
color_red <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[1]
color_green <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[3]
color_yellow <- RColorBrewer::brewer.pal(n = 7, name = "Set2")[6]
color_light <- RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")[3]
color_dark <- RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")[6]

file2write <- paste0(dir_out, "increased_proteins", ".png")
png(file2write, height = 1000, width = 800, res = 150)
p <- plot(plot_data_euler, labels = T, quantities = list(fontsize = 10), fills = c(color_light, color_dark, color_light))
grid.draw(p)
dev.off()

