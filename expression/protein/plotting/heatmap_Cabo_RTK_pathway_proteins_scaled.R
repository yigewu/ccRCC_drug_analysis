# Yige Wu @ WashU 2021 Jan

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

# input dependencies ------------------------------------------------------
## input the protein data
protein_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/preprocess/normalize_DIA_protein_data/20210111.v1/RCC_PDX.DIA_Protein.Log2.QuantileNormalized.20210111.v1.tsv")
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")

# make data matrix --------------------------------------------------------
## Cabozantinib: VEGFR2 inhibitor, also inhibits c-Met, Ret, Kit, Flt-1/3/4, Tie2, and AXL
met_related_genes <- c("HGF", "MET", "AXL")
vegfr_genes <- c("FLT1", "KDR", "FLT3", "FLT4", "NRP1", "NRP2")
vegf_genes <- c("VEGFA", "VEGFB", "VEGFC", "VEGFD", "VEGFE")
other_cabo_related_genes <- c("KIT", "RET", "NTRK2", "TEK")

## sum up the genes
genes_plot <- c(met_related_genes, vegfr_genes, vegf_genes, other_cabo_related_genes)
genes_plot <- unique(genes2plot)
genes_plot <- c(genes_plot, "Flt1", "Flt4")
genes_plot[!(genes_plot %in% protein_df$PG.Genes)]
## filter by gene
plot_data_df <- protein_df[protein_df$PG.Genes %in% genes_plot,]
plot_data_raw_mat <- as.matrix(plot_data_df[,5:ncol(plot_data_df)])
rownames(plot_data_raw_mat) <- plot_data_df$PG.Genes

## scale by row
plot_data_mat <- t(apply(plot_data_raw_mat, 1, scale))
rownames(plot_data_mat) <- plot_data_df$PG.Genes
colnames(plot_data_mat) <- colnames(plot_data_raw_mat)

# matrix row and column ids -----------------------------------------------
sampleids_mat <- colnames(plot_data_raw_mat)


# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color function for heatmap body colors
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
summary(as.vector(plot_data_mat))
colors_heatmapbody = circlize::colorRamp2(c(-1, 
                                            0, 
                                            1), 
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


# make row annotation -----------------------------------------------------
orig_avgexp_vec <- rowMeans(x = plot_data_raw_mat, na.rm = T)
row_anno_obj <- rowAnnotation(Unscaled_Expression = anno_simple(x = orig_avgexp_vec, col = colors_unscaledexp))

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
                                 

# plot  ------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = colors_heatmapbody,
                             na_col = color_na,
                             ## row args
                             left_annotation = row_anno_obj,
                             row_names_side = "left", show_row_dend = F,
                             ## column args
                             column_names_side = "top",
                             top_annotation = top_col_anno,
                             show_heatmap_legend = F)
p


# make legend -------------------------------------------------------------
annotation_lgd = list(
  Legend(labels = names(colors_treatment), 
         title = "Treatment", 
         legend_gp = gpar(fill = colors_treatment)),
  Legend(labels = names(colors_treatmentlength), 
         title = "Treatment length", 
         legend_gp = gpar(fill = colors_treatmentlength)))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "heatmap.png")
png(file2write, width = 1500, height = 1000, res = 150)
draw(object = p, 
     annotation_legend_side = "right", annotation_legend_list = annotation_lgd)
dev.off()
