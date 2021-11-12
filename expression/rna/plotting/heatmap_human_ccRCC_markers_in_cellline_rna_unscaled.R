# Yige Wu @ WashU 2021 Jan

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
source("./ccRCC_drug_analysis/plotting.R")
library(ComplexHeatmap)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the protein data
rna_exp_df <- fread("./Resources/Bulk_Processed_Data/Data_Files/batch8_cellline/geneExp/kallisto.tpm.gene_level.tsv", data.table = F)
genes_process_df <- fread(data.table = F, input = "../ccRCC_snRNA/Resources/Analysis_Results/findmarkers/tumor_specific_markers/overlap_tumor_vs_pt_DEGs_w_tumor_vs_other_DEGs/20210712.v1/ccRCC_markers.Surface.20210712.v1.tsv")

# set parameters ----------------------------------------------------------
## sum up the genes
genes2plot <- genes_process_df$Gene
genes2plot <- genes2plot[!(genes2plot %in% c("PIK3CB", "ARHGEF28", "PTGER3", "PARD3", "GNG12", "EFNA5", "SPIRE1", "LIFR", "PKP4", "SORBS1", "PTPRM", "FBXO16", "PAM"))]
genes2plot <- c("VEGFA", "GSTP1", genes2plot)
## set plotting metrics
num_nonna <- 0
row_fontsize <- 9

# make data matrix --------------------------------------------------------
## make the matrix to plot the heatmap
rna_exp_filtered_df <- rna_exp_df %>%
  filter(Name %in% genes2plot)
plot_data_raw_mat <- rna_exp_filtered_df %>%
  select(-Name)
plot_data_raw_mat <- as.matrix(plot_data_raw_mat)
## add row names
rownames(plot_data_raw_mat) <- rna_exp_filtered_df$Name
## transform the value
plot_data_raw_mat <- log2(plot_data_raw_mat+1)
plot_data_raw_mat <- plot_data_raw_mat[intersect(genes2plot, rownames(plot_data_raw_mat)),]
## scale by row
plot_data_mat <- t(apply(plot_data_raw_mat, 1, scale))
rownames(plot_data_mat) <- rownames(plot_data_raw_mat)
colnames(plot_data_mat) <- colnames(plot_data_raw_mat)

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color for treatment type
RColorBrewer::display.brewer.pal(name = "Set1", n = 8)
colors_treatment <- c("grey50", RColorBrewer::brewer.pal(n = 5, name = "Set1")[c(1,2,4,3)])
names(colors_treatment) <- c("Baseline", "Treated.Cabo", "Treated.Sap", "Treated.Cabo+Sap", "Control")
## make color for treatment length
colors_treatmentlength <- c(RColorBrewer::brewer.pal(name = "Set1", n = 8)[c(7, 5)], "grey50")
names(colors_treatmentlength) <- c("2 month", "1 month", "0 month")
## make colors for the original unscaled expression
colors_unscaledexp = circlize::colorRamp2(0:8, 
                                          RColorBrewer::brewer.pal(name = "RdPu", n = 9))
## make colors for species
colors_genespecies <- RColorBrewer::brewer.pal(name = "Dark2", n = 3)
names(colors_genespecies) <- c("Human", "Mouse", "Ambiguous")
## make colors for druggable
colors_yesno <- c("black", "grey80")
names(colors_yesno) <- c("TRUE", "FALSE")

# make row annotation -----------------------------------------------------
## get if each gene is druggable
isdruggable_vec <- mapvalues(x = rownames(plot_data_mat), from = genes_process_df$Gene, to = as.vector(genes_process_df$is_druggable))
isdruggable_vec[rownames(plot_data_mat) == "ENPP3"] <- "TRUE"
row_anno_obj <- rowAnnotation(is_druggable = anno_simple(x = isdruggable_vec, col = colors_yesno), 
                              annotation_name_side = "bottom")

## make color function for heatmap body colors
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
colors_heatmapbody <- colors_unscaledexp

# make column labels ------------------------------------------------------
colnames_plot <- colnames(plot_data_mat)
collabels_plot <- str_split_fixed(string = colnames_plot, pattern = "\\.", n = 2)[,1]

# make column split -------------------------------------------------------
col_split_vec <- ifelse(grepl(x = collabels_plot, pattern = "hk\\-2", ignore.case = T)| grepl(x = collabels_plot, pattern = "hek\\-293", ignore.case = T), "control", "ccRCC")

# make heatmap ------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_raw_mat, col = colors_heatmapbody,
                             ## row
                             cluster_rows = T, cluster_row_slices = F,
                             show_row_dend = F,
                             row_names_gp = grid::gpar(fontsize = row_fontsize), 
                             row_names_side = "left", 
                             right_annotation = row_anno_obj, 
                             row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 12),
                             ## column
                             column_split = col_split_vec, cluster_column_slices = T,
                             # column_title_rot = 90, column_title_gp = grid::gpar(fontsize = 10), top_annotation = top_col_anno,
                             # top_annotation = top_col_anno,
                             show_column_names = T, column_labels = collabels_plot,
                             cluster_columns = F, 
                             show_column_dend = F, show_heatmap_legend = F)
p

# make legend -------------------------------------------------------------
annotation_lgd = list(
  Legend(col_fun = colors_heatmapbody, 
         title = "Unscaled gene\nexpression", 
         title_gp = gpar(fontsize = 10),
         labels_gp = gpar(fontsize = 10)))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "heatmap.png")
png(file2write, width = 600, height = 1000, res = 150)
draw(p, annotation_legend_side = "right", annotation_legend_list = annotation_lgd)  #Show the heatmap
dev.off()
