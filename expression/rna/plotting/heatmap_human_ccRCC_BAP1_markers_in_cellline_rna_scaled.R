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

# set parameters ----------------------------------------------------------
## sum up the genes
genes2plot <- c("CA9", "CP", "MXI1", "KLF9", "SREBF2", "PCSK6", "PTPRJ", "CES3")
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
# ## filter columns
# plot_data_mat <- plot_data_mat[, c("dr-786-o-vhl-neg_rna.GATCTATC-ATGAGGCT",
#                                    "dr-hek-293t-rna.CATAATAC-TTCTAACG",
#                                    "dr-hk-2-rna.TGCGGCGT-CCTCGGTA",
#                                    "dr-rcc-4-vhl-pos_rna.AGCTCGCT-GCAGAATC",
#                                    "dr-rcc-4-vhl-neg_rna.CGGAACTG-CACTACGA")]

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color function for heatmap body colors
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
summary(as.vector(plot_data_mat))
colors_heatmapbody <- circlize::colorRamp2(c(quantile(plot_data_mat, 0.1, na.rm=T), 
                                             quantile(plot_data_mat, 0.5, na.rm=T), 
                                             quantile(plot_data_mat, 0.9, na.rm=T)),c(color_blue, "white", color_red))
colors_heatmapbody <- circlize::colorRamp2(c(-1.5, 
                                             0, 
                                             1.5),c(color_blue, "white", color_red))

## make color for treatment type
RColorBrewer::display.brewer.pal(name = "Set1", n = 8)
colors_treatment <- c("grey50", RColorBrewer::brewer.pal(n = 5, name = "Set1")[c(1,2,4,3)])
names(colors_treatment) <- c("Baseline", "Treated.Cabo", "Treated.Sap", "Treated.Cabo+Sap", "Control")
## make color for treatment length
colors_treatmentlength <- c(RColorBrewer::brewer.pal(name = "Set1", n = 8)[c(7, 5)], "grey50")
names(colors_treatmentlength) <- c("2 month", "1 month", "0 month")
## make colors for the original unscaled expression
colors_unscaledexp = circlize::colorRamp2(0:4, 
                                          RColorBrewer::brewer.pal(name = "RdPu", n = 5))
## make colors for species
colors_genespecies <- RColorBrewer::brewer.pal(name = "Dark2", n = 3)
names(colors_genespecies) <- c("Human", "Mouse", "Ambiguous")

# make row annotation -----------------------------------------------------
orig_avgexp_vec <- rowMeans(x = plot_data_raw_mat, na.rm = T)
row_anno_obj <- rowAnnotation(unscaled_expression = anno_simple(x = orig_avgexp_vec, col = colors_unscaledexp), 
                              annotation_name_side = "bottom")

# make column annotation --------------------------------------------------
treatment_vec <- mapvalues(x = sample_ids, from = data_status_df$SampleID.AcrossDataType, to = as.vector(data_status_df$ShortTag))
treatmentlength_vec <- mapvalues(x = sample_ids, from = data_status_df$SampleID.AcrossDataType, to = as.vector(data_status_df$Treatment.Month))
treatmentlength_vec <- paste0(treatmentlength_vec, " month")
## top column annotation object
top_col_anno = HeatmapAnnotation(TreatmentLength = anno_simple(x = treatmentlength_vec,
                                                               simple_anno_size = unit(3, "mm"),
                                                               col = colors_treatmentlength[treatmentlength_vec]),
                                 Cabozantinib = anno_simple(x = treatment_vec,
                                                            simple_anno_size = unit(3, "mm"),
                                                            col = colors_treatment), annotation_name_side = "left")

# make column labels ------------------------------------------------------
colnames_plot <- colnames(plot_data_mat)
collabels_plot <- str_split_fixed(string = colnames_plot, pattern = "\\.", n = 2)[,1]
collabels_plot <- gsub(x = collabels_plot, pattern = "_rna", replacement = "")
collabels_plot <- gsub(x = collabels_plot, pattern = "\\-rna", replacement = "")
collabels_plot <- gsub(x = collabels_plot, pattern = "dr\\-", replacement = "")

# make column split -------------------------------------------------------
col_split_vec <- ifelse(grepl(x = collabels_plot, pattern = "hk\\-2", ignore.case = T)| grepl(x = collabels_plot, pattern = "hek\\-293", ignore.case = T), "control", "ccRCC")

# make heatmap ------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, col = colors_heatmapbody,
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
         title = "Scaled gene\nexpression", 
         title_gp = gpar(fontsize = 10),
         labels_gp = gpar(fontsize = 10),
         # legend_width = unit(3, "cm"),
         legend_height = unit(1.5, "cm")),
  Legend(col_fun = colors_unscaledexp, 
         title = "Unscaled gene\nexpression", 
         title_gp = gpar(fontsize = 10),
         labels_gp = gpar(fontsize = 10),
         # legend_width = unit(3, "cm"),
         legend_height = unit(1.5, "cm")))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "heatmap.png")
png(file2write, width = 600, height = 800, res = 150)
draw(p, annotation_legend_side = "right", annotation_legend_list = annotation_lgd)  #Show the heatmap
dev.off()

file2write <- paste0(dir_out, "heatmap.pdf")
pdf(file2write, width = 3.5, height = 3, useDingbats = F)
draw(p, annotation_legend_side = "right", annotation_legend_list = annotation_lgd)  #Show the heatmap
dev.off()
