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
rna_exp_df <- fread("./Data_Freeze/v1.dataFreeze.washU_rcc/3.geneExp/v3.20210116/datafreeze.v3.kallisto.geneExp.protein_coding.tsv", data.table = F)
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")
## input the DEGs
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/examine_degs/filter_mouse_endothelialcells_RESL5_high/20210330.v1/Mouse.Endothelial.RESL5_vs_RESL10_CT.High.20210330.v1.tsv")
data_status_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/0.sample_info/v3.20210116/RCC_PDX_Samples.20210115.v2.xlsx")

# set parameters ----------------------------------------------------------
## sum up the genes
genes2plot <- unique(deg_df$deg_gene_symbol)
## set plotting metrics
num_nonna <- 0
row_fontsize <- 9

# get ids in order --------------------------------------------------------
## sort data status
data_status_filtered_df <- data_status_df %>%
  filter((Analysis_ID %in% colnames(rna_exp_df))) %>%
  filter(ShortTag %in% c("Baseline", "Control", "Treated.Cabo", "Treated.Cabo+Sap", "Treated.Sap"))
  arrange(ModelID, ShortTag)
analysis_ids <- data_status_filtered_df$Analysis_ID
analysis_ids
sample_ids <- mapvalues(x = analysis_ids, from = data_status_df$Analysis_ID, to = data_status_df$SampleID.AcrossDataType)
sample_ids
model_ids <- mapvalues(x = analysis_ids, from = data_status_df$Analysis_ID, to = data_status_df$ModelID)
model_ids
passages <- mapvalues(x = analysis_ids, from = data_status_df$Analysis_ID, to = data_status_df$NCI_Passage)
passages

# make data matrix --------------------------------------------------------
## make the matrix to plot the heatmap
rna_exp_filtered_df <- rna_exp_df %>%
  filter(Name %in% genes2plot)
plot_data_raw_mat <- rna_exp_filtered_df %>%
  select(-Name)
plot_data_raw_mat <- as.matrix(plot_data_raw_mat)
## add row names
rownames(plot_data_raw_mat) <- rna_exp_filtered_df$Name
## sort by column
plot_data_raw_mat <- plot_data_raw_mat[, analysis_ids]
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
## make color function for heatmap body colors
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
colors_heatmapbody <- circlize::colorRamp2(c(quantile(plot_data_mat, 0.1, na.rm=T), quantile(plot_data_mat, 0.5, na.rm=T), quantile(plot_data_mat, 0.9, na.rm=T)),c(color_blue, "white", color_red))

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

# make row annotation -----------------------------------------------------
orig_avgexp_vec <- rowMeans(x = plot_data_raw_mat, na.rm = T)
row_anno_obj <- rowAnnotation(Unscaled_Expression = anno_simple(x = orig_avgexp_vec, col = colors_unscaledexp), 
                              annotation_name_side = "top")


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



# make column split -------------------------------------------------------
col_split_vec <- model_ids
col_split_factor <- factor(x = col_split_vec, levels = c("RESL10", "RESL4", "RESL12", "RESL11", "RESL5", "RESL3", "RESL6", "RESL8", "RESL9"))
col_split_factor

# make heatmap ------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, col = colors_heatmapbody,
                             ## row
                             cluster_rows = T,
                             show_row_dend = F,
                             row_names_gp = grid::gpar(fontsize = row_fontsize), 
                             row_names_side = "left", right_annotation = row_anno_obj,
                             ## column
                             column_split = col_split_factor, 
                             column_title_rot = 90, column_title_gp = grid::gpar(fontsize = 10), top_annotation = top_col_anno,
                             # top_annotation = top_col_anno,
                             show_column_names = F, cluster_column_slices = T,
                             cluster_columns = T, 
                             show_column_dend = F, show_heatmap_legend = F)
p

# make legend -------------------------------------------------------------
annotation_lgd = list(
  Legend(labels = names(colors_treatment), 
         title = "Treatment", 
         legend_gp = gpar(fill = colors_treatment)),
  Legend(labels = names(colors_treatmentlength), 
         title = "Treatment length", 
         legend_gp = gpar(fill = colors_treatmentlength)),
  Legend(labels = names(colors_genespecies), 
         title = "Protein species", 
         legend_gp = gpar(fill = colors_genespecies)),
  Legend(col_fun = colors_heatmapbody, 
         title = "Scaled gene\nexpression", 
         title_gp = gpar(fontsize = 10),
         labels_gp = gpar(fontsize = 10),
         legend_width = unit(3, "cm"),
         legend_height = unit(3, "cm"),
         direction = "vertical"),
  Legend(col_fun = colors_unscaledexp, 
         title = "Unscaled gene\nexpression", 
         title_gp = gpar(fontsize = 10),
         labels_gp = gpar(fontsize = 10),
         legend_width = unit(3, "cm"),
         legend_height = unit(3, "cm"),
         direction = "vertical"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "heatmap.png")
png(file2write, width = 1800, height = 1200, res = 150)
p <- ComplexHeatmap::draw(p, annotation_legend_side = "right", annotation_legend_list = annotation_lgd)  #Show the heatmap
p
dev.off()
