# Yige Wu @WashU Dec 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
source("./ccRCC_drug_analysis/plotting.R")
library(readxl)
library(ComplexHeatmap)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies --------------------------------------------------
## input formatted summary data
summary_df <- read_excel(path = "./Resources/Treatment_Lists/Treatment_Summary/20200525_RCC single drug treatment summary_YW_simplified.xlsx", sheet = "all")

# make plot data ----------------------------------------------------------
plot_data_df <- summary_df[, c("Name", paste0("RESL", c(5, 10, 12, 4, 3, 11)))]
plot_data_df <- plot_data_df[rowSums(plot_data_df[,-1] != "Didn’t test") > 0,]
plot_data_mat <- as.matrix(x = plot_data_df[,-1])
rownames(plot_data_mat) <- plot_data_df$Name
plot_data_mat[plot_data_mat == "Less effective"] <- "Effective"
plot_data_mat[plot_data_mat == "Didn’t test"] <- "Didn't test"

plot_data_mat <- plot_data_mat[c("Cabozantinib",
                                 "Sapanisertib",
                                 "Panobinostat",
                                 "Sunitinib",
                                 "Abemaciclib",
                                 "Selumetinib",
                                 "Birinapant",
                                 "Beta-hydroxybutyrate",
                                 "Acriflavine hydrochloride",
                                 "Losartan"),]

# make colors -------------------------------------------------------------
colors_heatmapbody <- structure(c("black", "grey50", "red", "orange"), names = c("Didn’t test", "Not effective", "Effective", "Less effective"))
colors_heatmapbody <- structure(c("black", "grey50", "red"), names = c("Didn't test", "Not effective", "Effective"))

# plot --------------------------------------------------------------------
fontsize_plot <- 15
p <- Heatmap(matrix = plot_data_mat, col = colors_heatmapbody, show_heatmap_legend = F, 
             row_names_gp = gpar(fontsize = fontsize_plot), column_names_gp = gpar(fontsize = fontsize_plot), cluster_rows = F)
annotation_lgd = list(
  Legend(labels = names(colors_heatmapbody), labels_gp = gpar(fontsize = fontsize_plot),
         title = "Response", title_gp = gpar(fontsize = fontsize_plot),
         legend_gp = gpar(fill = colors_heatmapbody)))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "heatmap.pdf")
pdf(file2write, width = 5, height = 3.5, useDingbats = F)
draw(p, annotation_legend_side = "right", annotation_legend_list = annotation_lgd)
dev.off()
