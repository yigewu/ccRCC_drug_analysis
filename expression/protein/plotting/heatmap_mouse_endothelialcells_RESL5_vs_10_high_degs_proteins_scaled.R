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
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the protein data
protein_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/preprocess/impute_protein_data_using_DreamAI/20210201.v1/RCC_PDX.DIA_Protein.Log2.QuantileNormalized.RegImpute.tsv")
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")
## input the anova results
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/examine_degs/filter_mouse_endothelialcells_RESL5_high/20210330.v1/Mouse.Endothelial.RESL5_vs_RESL10_CT.High.20210330.v1.tsv")

# set parameters ----------------------------------------------------------
## sum up the genes
proteins_plot <- deg_df$PG.ProteinNames[!is.na(deg_df$PG.ProteinNames)]
proteins_plot <- unique(proteins_plot)
length(proteins_plot)

# make data matrix --------------------------------------------------------
## filter by gene
plot_data_df <- protein_df[protein_df$PG.ProteinNames %in% proteins_plot,]
plot_data_raw_mat <- as.matrix(plot_data_df[,5:ncol(plot_data_df)])
rownames(plot_data_raw_mat) <- plot_data_df$PG.ProteinNames

## scale by row
plot_data_mat <- t(apply(plot_data_raw_mat, 1, scale))
rownames(plot_data_mat) <- plot_data_df$PG.ProteinNames
colnames(plot_data_mat) <- colnames(plot_data_raw_mat)

# matrix row and column ids -----------------------------------------------
sampleids_mat <- colnames(plot_data_mat)
proteinnames_mat <- rownames(plot_data_mat)

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color function for heatmap body colors
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
summary(as.vector(plot_data_mat))
# colors_heatmapbody = circlize::colorRamp2(c(rev(seq(-2, -5, -1)), rev(seq(-0.2, -1, -0.2)), 
#                                   0,
#                                   seq(0.2, 1, 0.2), seq(2, 5, 1)),
#                                 c(rev(brewer.pal(n = 9, name = "Blues")), "white", brewer.pal(n = 9, name = "YlOrRd")))
colors_heatmapbody = circlize::colorRamp2(c(-1.25,
                                            0,
                                            1.25),
                                          c(color_blue, "white", color_red))
## make color for treatment type
RColorBrewer::display.brewer.pal(name = "Set1", n = 8)
colors_treatment <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[c(1,2,4,3)]
names(colors_treatment) <- c("Cabo", "Sap", "Cabo+ Sap", "Con")
## make color for treatment length
colors_treatmentlength <- RColorBrewer::brewer.pal(name = "Set1", n = 8)[c(7, 5)]
names(colors_treatmentlength) <- c("2 month", "1 month")
## make colors for the original unscaled expression
colors_unscaledexp = circlize::colorRamp2(14:22, 
                                          RColorBrewer::brewer.pal(name = "RdPu", n = 9))
## make colors for species
colors_genespecies <- RColorBrewer::brewer.pal(name = "Dark2", n = 3)
names(colors_genespecies) <- c("Human", "Mouse", "Ambiguous")

# make row annotation -----------------------------------------------------
orig_avgexp_vec <- rowMeans(x = plot_data_raw_mat, na.rm = T)
## map gene to human/mouse
ishuman_vec <- grepl(pattern = "HUMAN", x = proteinnames_mat)
ismouse_vec <- grepl(pattern = "MOUSE", x = proteinnames_mat)
protein_species_vec <- ifelse(ishuman_vec & !ismouse_vec, "Human",
                           ifelse(ismouse_vec & !ishuman_vec, "Mouse", "Ambiguous"))
row_anno_obj <- rowAnnotation(Unscaled_Expression = anno_simple(x = orig_avgexp_vec, col = colors_unscaledexp), 
                              Species = anno_simple(x = protein_species_vec, col = colors_genespecies),
                              annotation_name_side = "top")

# make column annotation --------------------------------------------------
treatment_vec <- mapvalues(x = sampleids_mat, from = sampleinfo_df$`Sample ID`, to = as.vector(sampleinfo_df$Treatment))
treatmentlength_vec <- mapvalues(x = sampleids_mat, from = sampleinfo_df$`Sample ID`, to = as.vector(sampleinfo_df$Treatment_length))

top_col_anno = HeatmapAnnotation(TreatmentLength = anno_simple(x = treatmentlength_vec,
                                                               simple_anno_size = unit(3, "mm"),
                                                               # gp = gpar(col = "black"),
                                                               col = colors_treatmentlength[treatmentlength_vec]),
                                 Treatment = anno_simple(x = treatment_vec,
                                                         simple_anno_size = unit(3, "mm"),
                                                         # gp = gpar(col = "black"),
                                                         col = colors_treatment[treatment_vec]),
                                 annotation_name_side = "left")

# make column split -------------------------------------------------------
col_split_vec <- str_split_fixed(string = sampleids_mat, pattern = "_", 3)[,1]
col_split_factor <- factor(x = col_split_vec, levels = c("RESL5", "RESL4", "RESL10", "RESL12", "RESL3", "RESL11", "RESL6", "RESL8", "RESL9"))
col_split_factor


# plot  ------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = colors_heatmapbody,
                             na_col = color_na,
                             ## row args
                             right_annotation = row_anno_obj,
                             show_row_names = T,
                             row_names_side = "left",
                             show_row_dend = F, 
                             # row_km = 6, row_km_repeats = 100,
                             ## column args
                             column_split = col_split_factor, show_column_names = F,
                             # column_names_side = "top", 
                             top_annotation = top_col_anno, 
                             cluster_column_slices = F, cluster_columns = F,
                             show_column_dend = F, column_dend_height = unit(2, "cm"), 
                             # column_km = 6, column_km_repeats = 200,
                             show_heatmap_legend = F)

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
         title = "Scaled protein\nabundance", 
         title_gp = gpar(fontsize = 10),
         labels_gp = gpar(fontsize = 10),
         legend_width = unit(3, "cm"),
         legend_height = unit(3, "cm"),
         direction = "vertical"),
  Legend(col_fun = colors_unscaledexp, 
         title = "Unscaled protein\nabundance", 
         title_gp = gpar(fontsize = 10),
         labels_gp = gpar(fontsize = 10),
         legend_width = unit(3, "cm"),
         legend_height = unit(3, "cm"),
         direction = "vertical"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "heatmap.png")
png(file2write, width = 1800, height = 2000, res = 150)
p <- ComplexHeatmap::draw(p, annotation_legend_side = "right", annotation_legend_list = annotation_lgd)  #Show the heatmap
print(p)
dev.off()
