# Yige Wu @ WashU 2021 Jan

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
packages = c(
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2",
  "ComplexHeatmap"
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
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the protein data
exp_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/unite_protein_phospho_diff_treated_vs_control_bylog2intensity/20220110.v1/Protein_Phospho_Diff_Log2Intensity.Treated_vs_Control.20220110.v1.tsv")
## input the pathway score
pathwayscores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/calculate_pathway_score_treatedvscontrol_bydiff_log2intensity/20220112.v1/Pathway_scores_bysample.20220112.v1.tsv")
pathwayscores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/calculate_pathway_score_treatedvscontrol_bydiff_log2intensity/20220112.v3/Pathway_scores_bysample.20220112.v3.tsv")

## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")

# get proteins to plot ----------------------------------------------------------
atrophygenes_plot = c("Foxo1", "Socs3", "Stat3", "Acvr2b", "Fbxo32", "Trim63")
immunegenes_mouse_plot = c("Cd40lg", "Egf", "Ccl11", "Fgf2", "Flt3lg", "Cx3cl1", "Csf3", "Csf2", "Cxcl1", "Cxcl2", "Cxcl3", 
                           "Ifna2", "Ifng", "Il1a", "Il1b", "Il1rn", "Il2", "Il3", "Il4", "Il5", "Il6", "Il7", "Cxcl8", 
                           "Il9", "Il10","Il12b","Il12a","Il12b","Il13","Il15","Il17a","Cxcl10","Ccl2","Ccl7","Ccl22","Ccl3","Ccl4","Tgfa","Tnf","Lta","Vegfa")
immunegenes_human_plot <- c("CD40LG", "EGF", "CCL11", "FGF2", "FLT3LG", "CX3CL1", "CSF3", "CSF2", "CXCL1", "CXCL2", "CXCL3", 
                            "IFNA2", "IFNG", "IL1A", "IL1B", "IL1RN", "IL2", "IL3", "IL4", "IL5", "IL6", "IL7", "CXCL8", 
                            "IL9", "IL10","IL12B","IL12A","IL12B","IL13","IL15","IL17A","CXCL10","CCL2","CCL7","CCL22","CCL3","CCL4","TGFA","TNF","LTA","VEGFA")
proteins_plot <- exp_df$PG.Gene[exp_df$PG.Gene %in% c(atrophygenes_plot, 
                                                      immunegenes_mouse_plot, immunegenes_human_plot)]
proteins_plot = unique(proteins_plot)
length(proteins_plot)
exp_df$ID <- paste0(exp_df$PG.Gene, "_", exp_df$PTM_Name)

# for (treatment_group in c("Cabo")) {
for (treatment_group in c("Cabo", "Sap", "Cabo+ Sap")) {
  # make data matrix --------------------------------------------------------
  ## filter by gene
  plot_data_df <- exp_df %>%
    filter(PG.Gene %in% proteins_plot) %>%
    filter(PTM_Name == "Protein")
  sampleinfo_df <- sampleinfo_df %>%
    mutate(model_id = str_split_fixed(string = `Sample ID`, pattern = "_", n = 3)[,1]) %>%
    mutate(sample_id = paste0(model_id, "_", Treatment, "_", gsub(x = Treatment_length, pattern = " ", replacement = "")))
  sampleinfo_plot_df <- sampleinfo_df %>%
    filter(Treatment == treatment_group) %>%
    filter(Treatment_length == "1 month") %>%
    arrange(Treatment, Treatment_length, model_id)
  # sample_ids <- colnames(plot_data_df); sample_ids <- sample_ids[grepl(pattern = "RESL", x = sample_ids)]; #sample_ids <- sample_ids[grepl(pattern = "_Cabo_", x = sample_ids)]
  sample_ids <- sampleinfo_plot_df$sample_id
  plot_data_mat <- as.matrix(plot_data_df[,sample_ids])
  rownames(plot_data_mat) <- plot_data_df$ID
  
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
  # RColorBrewer::display.brewer.pal(name = "Set1", n = 8)
  colors_treatment <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[c(1,2,4,3)]
  names(colors_treatment) <- c("Cabo", "Sap", "Cabo+ Sap", "Con")
  ## make color for treatment length
  colors_treatmentlength <- RColorBrewer::brewer.pal(name = "Set1", n = 8)[c(7, 5)]
  names(colors_treatmentlength) <- c("2month", "1month")
  colors_bymodel <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[-6]
  names(colors_bymodel) <- c("RESL5", "RESL10", "RESL12", "RESL4", "RESL11", "RESL3")
  colors_pathwayscore = circlize::colorRamp2(seq(-0.5, 0.5, 0.1),
                                             rev(brewer.pal(n = 11, name = "BrBG")))
  
  # # make row annotation -----------------------------------------------------
  # row_anno_obj <- rowAnnotation(Unscaled_Expression = anno_simple(x = orig_avgexp_vec, col = colors_unscaledexp),
  #                               Species = anno_simple(x = protein_species_vec, col = colors_genespecies),
  #                               annotation_name_side = "top")
  
  # make row split --------------------------------------------------
  # row_split_vec <- mapvalues(x = proteinnames_mat, from = pathway2members_df$ID, to = as.vector(pathway2members_df$Type_of_regulation))
  
  # make column annotation --------------------------------------------------
  treatment_vec <- str_split_fixed(string = sampleids_mat, pattern = "_", n = 3)[,2]
  treatmentlength_vec <- str_split_fixed(string = sampleids_mat, pattern = "_", n = 3)[,3]
  modelid_vec <- str_split_fixed(string = sampleids_mat, pattern = "_", n = 3)[,1]
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
                                   annotation_name_side = "left")
  
  # make column split -------------------------------------------------------
  col_split_vec <- str_split_fixed(string = sampleids_mat, pattern = "_", n = 3)[,2]
  
  # plot  ------------------------------------------------------------
  p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                               col = colors_heatmapbody,
                               na_col = color_na,
                               ## row args
                               # right_annotation = row_anno_obj,
                               # row_split = row_split_vec,
                               show_row_names = T,
                               row_names_side = "left",
                               show_row_dend = F, cluster_rows = F,
                               # row_km = 6, row_km_repeats = 100,
                               ## column args
                               column_names_side = "top", 
                               # column_split = treatment_vec, 
                               column_split = modelid_vec, 
                               cluster_columns = F, show_column_names = F,
                               # top_annotation = top_col_anno,
                               show_column_dend = F, #column_dend_height = unit(2, "cm"),
                               show_heatmap_legend = F)
  p
  # make legend -------------------------------------------------------------
  annotation_lgd = list(
    Legend(col_fun = colors_heatmapbody, 
           title = "Protein abundance difference (treated vs. control)", 
           title_gp = gpar(fontsize = 10),
           labels_gp = gpar(fontsize = 10),
           legend_width = unit(3, "cm"),
           legend_height = unit(3, "cm"),
           direction = "horizontal"))
  
  # write output ------------------------------------------------------------
  file2write <- paste0(dir_out, "heatmap.", treatment_group, ".png")
  png(file2write, width = 850, height = 450, res = 150)
  draw(p, annotation_legend_side = "bottom", annotation_legend_list = annotation_lgd)  #Show the heatmap
  dev.off()
  
}
