# Yige Wu @ WashU 2022 Jan

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
exp_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/unite_protein_phospho_diff_treated_vs_control_bylog2intensity/20220110.v1/Protein_Phospho_Diff_Log2Intensity.Treated_vs_Control.20220110.v1.tsv")
## input the pathway score
# pathwayscores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/calculate_pathway_score_treatedvscontrol_bydiff_log2intensity/20220112.v1/Pathway_scores_bysample.20220112.v1.tsv")
pathwayscores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/calculate_pathway_score_treatedvscontrol_bydiff_log2intensity/20220112.v3/Pathway_scores_bysample.20220112.v3.tsv")

## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")

# get proteins to plot ----------------------------------------------------------
pathwayname_plot <- "PI3K_Akt"
pathway2members_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Gene_Lists/Pathway_score_members.011222.xlsx", sheet = pathwayname_plot)
pathway2members_df <- pathway2members_df %>%
  mutate(ID = paste0(Gene_symbol, "_", ifelse(Is_phospho == "Yes", Site_phospho, "Protein")))
# proteins_plot <- paste0(proteins_plot, "_Protein")
proteins_plot <- pathway2members_df$ID
length(proteins_plot)
exp_df$ID <- paste0(exp_df$PG.Gene, "_", exp_df$PTM_Name)

# make data matrix --------------------------------------------------------
## filter by gene
plot_data_df <- exp_df %>%
  filter(ID %in% proteins_plot)
sampleinfo_df <- sampleinfo_df %>%
  mutate(model_id = str_split_fixed(string = `Sample ID`, pattern = "_", n = 3)[,1]) %>%
  mutate(sample_id = paste0(model_id, "_", Treatment, "_", gsub(x = Treatment_length, pattern = " ", replacement = "")))
sampleinfo_plot_df <- sampleinfo_df %>%
  filter(Treatment != "Con") %>%
  arrange(Treatment, Treatment_length, model_id)
sample_ids <- colnames(plot_data_df); sample_ids <- sample_ids[grepl(pattern = "RESL", x = sample_ids)]; #sample_ids <- sample_ids[grepl(pattern = "_Cabo_", x = sample_ids)]
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
colors_heatmapbody = circlize::colorRamp2(seq(-1, 1, 0.2),
                                          rev(brewer.pal(n = 11, name = "RdBu")))
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

# # make row annotation -----------------------------------------------------
# row_anno_obj <- rowAnnotation(Unscaled_Expression = anno_simple(x = orig_avgexp_vec, col = colors_unscaledexp),
#                               Species = anno_simple(x = protein_species_vec, col = colors_genespecies),
#                               annotation_name_side = "top")

# make row split --------------------------------------------------
row_split_vec <- mapvalues(x = proteinnames_mat, from = pathway2members_df$ID, to = as.vector(pathway2members_df$Type_of_regulation))

# make column annotation --------------------------------------------------
treatment_vec <- str_split_fixed(string = sampleids_mat, pattern = "_", n = 3)[,2]
treatmentlength_vec <- str_split_fixed(string = sampleids_mat, pattern = "_", n = 3)[,3]
modelid_vec <- str_split_fixed(string = sampleids_mat, pattern = "_", n = 3)[,1]
pathwayscore_vec <- pathwayscores_df[pathwayscores_df$Pathway_name == pathwayname_plot, sampleids_mat];
pathwayscore_vec <- unlist(pathwayscore_vec)
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
                                 PathwayScore = anno_simple(x = pathwayscore_vec, col = colors_pathwayscore),
                                 annotation_name_side = "left")

# make column split -------------------------------------------------------
col_split_vec <- str_split_fixed(string = sampleids_mat, pattern = "_", n = 3)[,2]

# plot  ------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = colors_heatmapbody,
                             na_col = color_na,
                             ## row args
                             # right_annotation = row_anno_obj,
                             row_split = row_split_vec,
                             show_row_names = T,
                             row_names_side = "left",
                             show_row_dend = F, 
                             # row_km = 6, row_km_repeats = 100,
                             ## column args
                             column_names_side = "top", 
                             # column_split = treatment_vec, 
                             column_split = modelid_vec, 
                             cluster_columns = F, show_column_names = F,
                             top_annotation = top_col_anno,
                             show_column_dend = F, #column_dend_height = unit(2, "cm"),
                             show_heatmap_legend = F)
p
# make legend -------------------------------------------------------------
annotation_lgd = list(
  Legend(labels = names(colors_treatment),
         title = "Treatment",
         legend_gp = gpar(fill = colors_treatment)),
  Legend(labels = names(colors_treatmentlength), 
         title = "Treatment length", 
         legend_gp = gpar(fill = colors_treatmentlength)),
  # Legend(labels = names(colors_bymodel), 
  #        title = "Model", 
  #        legend_gp = gpar(fill = colors_bymodel)),
  Legend(col_fun = colors_pathwayscore, 
         title = "Pathway score difference", 
         title_gp = gpar(fontsize = 10),
         labels_gp = gpar(fontsize = 10),
         legend_width = unit(3, "cm"),
         legend_height = unit(3, "cm"),
         direction = "horizontal"),
  Legend(col_fun = colors_heatmapbody, 
         title = "Protein/phosphorylation\nabundance difference\n(treated vs. control)", 
         title_gp = gpar(fontsize = 10),
         labels_gp = gpar(fontsize = 10),
         legend_width = unit(3, "cm"),
         legend_height = unit(3, "cm"),
         direction = "horizontal"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "heatmap.png")
png(file2write, width = 1000, height = 400, res = 150)
draw(p, annotation_legend_side = "bottom", annotation_legend_list = annotation_lgd)  #Show the heatmap
dev.off()
# file2write <- paste0(dir_out, "heatmap.RDS")
# saveRDS(object = p, file = file2write, compress = T)
