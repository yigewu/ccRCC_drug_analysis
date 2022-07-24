# Yige Wu @ WashU 2022 Jan

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
  "ggplot2",
  "RColorBrewer",
  "ComplexHeatmap"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_drug_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the protein data
exp_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/unite_protein_phospho_bylog2intensity/20220120.v1/Protein_Phospho_Log2Intensity.20220120.v1.tsv")
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")
## input the pathway score
# pathwayscores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/calculate_pathway_score_treatedvscontrol_bydiff_log2intensity/20220112.v1/Pathway_scores_bysample.20220112.v1.tsv")
# pathwayscores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/calculate_pathway_score_treatedvscontrol_bydiff_log2intensity/20220112.v3/Pathway_scores_bysample.20220112.v3.tsv")

# get proteins to plot ----------------------------------------------------------
proteins_plot <- paste0(c(#"IGF2BP3", "PYCR1", "ERO1B", "DNAJC7", "CHMP6", "MRM3", "CLIC6", "COBLL1"#,
                          "TGFBI", "LRP1", "COL7A1", "TAGLN", "SCG2", "THBS1", "TFPI2"
                          ), "_Protein")
length(proteins_plot)

# make data matrix --------------------------------------------------------
## filter by gene
plot_data_df <- exp_df %>%
  filter(ID %in% proteins_plot)
sampleinfo_df <- sampleinfo_df %>%
  mutate(model_id = str_split_fixed(string = `Sample ID`, pattern = "_", n = 3)[,1]) %>%
  mutate(sample_id = paste0(model_id, "_", Treatment, "_", gsub(x = Treatment_length, pattern = " ", replacement = "")))
sampleinfo_plot_df <- sampleinfo_df %>%
  arrange(model_id, Treatment_length, Treatment)
sample_ids <- sampleinfo_plot_df$`Sample ID`
plot_data_raw_mat <- as.matrix(plot_data_df[,sample_ids])
## scale by row
plot_data_mat <- t(apply(plot_data_raw_mat, 1, scale))
rownames(plot_data_mat) <- plot_data_df$ID
colnames(plot_data_mat) <- sampleinfo_plot_df$sample_id

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
colors_heatmapbody = circlize::colorRamp2(seq(-2, 2, 0.4),
                                          rev(brewer.pal(n = 11, name = "RdBu")))
# colors_heatmapbody = circlize::colorRamp2(seq(-1, 1, 0.2),
#                                           rev(brewer.pal(n = 11, name = "RdBu")))
## make color for treatment type
RColorBrewer::display.brewer.pal(name = "Set1", n = 8)
colors_treatment <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[c(1,2,4,3)]
names(colors_treatment) <- c("Cabo", "Sap", "Cabo+ Sap", "Con")
## make color for treatment length
colors_treatmentlength <- RColorBrewer::brewer.pal(name = "Set1", n = 8)[c(7, 5)]
names(colors_treatmentlength) <- c("2month", "1month")
colors_bymodel <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[-6]
names(colors_bymodel) <- c("RESL5", "RESL10", "RESL12", "RESL4", "RESL11", "RESL3")
colors_pathwayscore = circlize::colorRamp2(seq(-0.5, 0.5, 0.1),
                                          rev(brewer.pal(n = 11, name = "BrBG")))
## make colors for the original unscaled expression
colors_unscaledexp = circlize::colorRamp2(14:18, 
                                          RColorBrewer::brewer.pal(name = "Greens", n = 5))
colors_typeofreg <- RColorBrewer::brewer.pal(n = 5, name = "Dark2")[3:4]
names(colors_typeofreg) <- c("Negative", "Positive")

# # make row annotation -----------------------------------------------------
# typeofregulation_vec <- mapvalues(x = proteinnames_mat, from = pathway2members_df$ID, to = as.vector(pathway2members_df$Type_of_regulation))
orig_avgexp_vec <- rowMeans(x = plot_data_raw_mat, na.rm = T)
row_anno_obj <- rowAnnotation(#TypeofRegulation = anno_simple(x = typeofregulation_vec, col = colors_typeofreg[typeofregulation_vec]),
                              Unscaled_Expression = anno_simple(x = orig_avgexp_vec, col = colors_unscaledexp), 
                              annotation_name_side = "top")

# make row split --------------------------------------------------
# row_split_vec <- mapvalues(x = proteinnames_mat, from = pathway2members_df$ID, to = as.vector(pathway2members_df$Pathway_name))
# # row_split_vec2 <- mapvalues(x = row_split_vec, from = )
# row_split_factor <- factor(x = row_split_vec, levels = pathwaynames_ordered)

# make column annotation --------------------------------------------------
treatment_vec <- str_split_fixed(string = sampleids_mat, pattern = "_", n = 3)[,2]
treatmentlength_vec <- str_split_fixed(string = sampleids_mat, pattern = "_", n = 3)[,3]
modelid_vec <- str_split_fixed(string = sampleids_mat, pattern = "_", n = 3)[,1]
# pathwayscore_vec <- pathwayscores_df[pathwayscores_df$Pathway_name == pathwayname_plot, sampleids_mat];
# pathwayscore_vec <- unlist(pathwayscore_vec)
top_col_anno = HeatmapAnnotation(TreatmentLength = anno_simple(x = treatmentlength_vec,
                                                               simple_anno_size = unit(3, "mm"),
                                                               # gp = gpar(col = "black"),
                                                               col = colors_treatmentlength[treatmentlength_vec]),
                                 Treatment = anno_simple(x = treatment_vec,
                                                         simple_anno_size = unit(3, "mm"),
                                                         # gp = gpar(col = "black"),
                                                         col = colors_treatment[treatment_vec]),
                                 # Model = anno_simple(x = modelid_vec,
                                 #                         simple_anno_size = unit(3, "mm"),
                                 #                         # gp = gpar(col = "black"),
                                 #                         col = colors_bymodel[modelid_vec]),
                                 # PathwayScore = anno_simple(x = pathwayscore_vec, col = colors_pathwayscore),
                                 annotation_name_side = "left")

# make column split -------------------------------------------------------
col_split_vec <- modelid_vec
# col_split_factor <- factor(x = col_split_vec, levels = c("RESL10", "RESL4", "RESL12", "RESL11", "RESL5", "RESL3", "RESL6", "RESL8", "RESL9"))
col_split_factor <- factor(x = col_split_vec, levels = c("RESL5", "RESL4", "RESL3", "RESL11", "RESL12", "RESL10"))
col_split_factor

# plot  ------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = colors_heatmapbody,
                             na_col = color_na,
                             ## row args
                             right_annotation = row_anno_obj,
                             # row_split = row_split_factor, 
                             # cluster_row_slices = F, row_title_rot = 0,
                             show_row_names = T,
                             row_names_side = "left",
                             show_row_dend = F, 
                             # row_km = 6, row_km_repeats = 100,
                             ## column args
                             column_names_side = "top", 
                             # column_split = treatment_vec, 
                             column_split = col_split_factor, column_title_rot = 90, cluster_column_slices = F,
                             cluster_columns = F, show_column_names = F,
                             top_annotation = top_col_anno,
                             show_column_dend = F, #column_dend_height = unit(2, "cm"),
                             show_heatmap_legend = F)
p
# make legend -------------------------------------------------------------
annotation_lgd = list(
  Legend(labels = names(colors_treatmentlength), 
         title = "Vehicle treatment\nlength", 
         legend_gp = gpar(fill = colors_treatmentlength)),
  Legend(labels = names(colors_treatment), 
         title = "Treatment type", 
         legend_gp = gpar(fill = colors_treatment)),
  # Legend(labels = names(colors_typeofreg),
  #        title = "Type of regulation\non the pathways",
  #        legend_gp = gpar(fill = colors_typeofreg)),
  # Legend(labels = names(colors_bymodel),
  #        title = "Model",
  #        legend_gp = gpar(fill = colors_bymodel)),
  Legend(col_fun = colors_heatmapbody, 
         title = "Protein/phosphorylation\nabundance", 
         title_gp = gpar(fontsize = 10),
         labels_gp = gpar(fontsize = 10),
         legend_width = unit(3, "cm"),
         legend_height = unit(3, "cm"),
         direction = "horizontal"),
  Legend(col_fun = colors_unscaledexp, 
         title = "Unscaled protein\nabundance", 
         title_gp = gpar(fontsize = 10),
         labels_gp = gpar(fontsize = 10),
         # legend_width = unit(3, "cm"),
         legend_height = unit(1.5, "cm"), 
         direction = "horizontal"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "heatmap.png")
png(file2write, width = 700, height = 900, res = 150)
draw(p, annotation_legend_side = "bottom", annotation_legend_list = annotation_lgd)  #Show the heatmap
dev.off()
file2write <- paste0(dir_out, "heatmap.pdf")
pdf(file2write, width = 6, height = 5, useDingbats = F)
draw(p, annotation_legend_side = "bottom", annotation_legend_list = annotation_lgd)  #Show the heatmap
dev.off()
