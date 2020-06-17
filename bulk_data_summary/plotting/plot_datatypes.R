# Yige Wu @WashU March 2020
## running on local
## for plotting the aliquot-pairwise correlation coefficients for averaged expression of all genes

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
## input all ccRCC samples data info
data_status_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_status/write_sn_bulk_data_status/20200616.v1/RCC_PDX_Related_Samples.PDXNet_B1_9.HTAN_B1.Bulk_SingleNuclei_Data_Status.20200616.v1.tsv")

# get ids in order --------------------------------------------------------
## sort data status
data_status_filtered_df <- data_status_df %>%
  filter(Group == "PDX") %>%
  # filter(!is.na(Analysis_ID.RNA) | !is.na(Analysis_ID.WES)) %>%
  arrange(ModelID, TreatmentLengthGroup, snRNA, Treated.Cab, Treated.Sap, Treated.ACF, Treated.Entinostat)
sample_ids <- data_status_filtered_df$SampleID.AcrossDataType
sample_ids
model_ids <- data_status_filtered_df$ModelID
model_ids
passages <- data_status_filtered_df$NCI_Passage
passages

# make data matrix for heatmap body ---------------------------------------
## make empty data frmae
plot_data_mat <- matrix(nc = length(sample_ids), nr = 0)
plot_data_mat %>% head()

# make top column annotation --------------------------------------------------
## filter by columns
top_col_anno_df <- data_status_filtered_df
## add row names
rownames(top_col_anno_df) <- sample_ids
## make color for sequencing status
RColorBrewer::display.brewer.pal(name = "Spectral", n = 11)
colors_seq_status <- c(RColorBrewer::brewer.pal(n = 11, name = "Spectral")[c(9,7,5)], "grey80")
names(colors_seq_status) <- c("Data Processed", "Data Processing", "Sequencing", "Not Sequenced")
## get unique color for different passages
passages <- top_col_anno_df$NCI_Passage
uniq_passages <- unique(passages)
uniq_passages
colors_uniq_passages <- RColorBrewer::brewer.pal(n = 5, name = "PuRd")[3:5]
names(colors_uniq_passages) <- c("P3", "P4", "P5")
## make color for  CNV
colors_cn <- rev(RColorBrewer::brewer.pal(n = 9, name = "PuOr"))[3:9]
names(colors_cn) <- 0:6
## make color for treatment length
treatmentlength <- top_col_anno_df$TreatmentLengthGroup
colors_treatmentlength <- RColorBrewer::brewer.pal(n = 4, name = "BuPu")[3:4]
names(colors_treatmentlength) <- c("2_Months", "1_Month")
## make color for treatment info
RColorBrewer::display.brewer.pal(name = "Set1", n = 5)
colors_treatment <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[c(1,3,2,4,5)]
names(colors_treatment) <- c("cab", "sap", "acf", "ent", "sun")
## get if it's control
is_control <- as.numeric(top_col_anno_df$TumorTissue == "PDX tumor;control")
colors_is_control <- c("1" = "grey", "0" = "white")
## get if it's cabo treated
is_cab_treated <- as.numeric(top_col_anno_df$Treated.Cab)
colors_is_cab_treated <- c(colors_treatment["cab"], "white")
names(colors_is_cab_treated) <- c(1,0)
## get if it's Sap treated
is_sap_treated <- as.numeric(top_col_anno_df$Treated.Sap)
colors_is_sap_treated <- c(colors_treatment["sap"], "white")
names(colors_is_sap_treated) <- c(1,0)
## get if it's ACF treated
is_acf_treated <- as.numeric(top_col_anno_df$Treated.ACF)
colors_is_acf_treated <- c(colors_treatment["acf"], "white")
names(colors_is_acf_treated) <- c(1,0)
## get if it's entinostat treated
is_ent_treated <- as.numeric(top_col_anno_df$Treated.Entinostat)
colors_is_ent_treated <- c(colors_treatment["ent"], "white")
names(colors_is_ent_treated) <- c(1,0)
## top column annotation object
top_col_anno = HeatmapAnnotation(TreatmentLength = anno_text(x = treatmentlength, 
                                                             location = 0.5, just = "center",
                                                             gp = gpar(fill = colors_treatmentlength[treatmentlength], col = "white", border = "black", fontsize = 8),
                                                             width = max_text_width(treatmentlength)*1.4),
                                 # ModelID = anno_text(x = model_ids, 
                                 #                     location = 0.5, just = "center",
                                 #                     gp = gpar(fill = uniq_model_colors[model_ids], col = "white", border = "black"),
                                 #                     width = max_text_width(model_ids)*1.2),
                                 # NCI_Passage = anno_text(x = passages, 
                                 #                         location = 0.5, just = "center",
                                 #                         gp = gpar(fill = colors_uniq_passages[passages], col = "white", border = "black"),
                                 #                         width = max_text_width(model_ids)*1.2),
                                 # CN.VHL.3p = anno_simple(x = top_col_anno_df$CN.VHL,
                                 #                         simple_anno_size = unit(3, "mm"), 
                                 #                         col = colors_cn),
                                 # CN.HIF1A.14q = anno_simple(x = top_col_anno_df$CN.HIF1A,
                                 #                            simple_anno_size = unit(3, "mm"), 
                                 #                            col = colors_cn),
                                 # CN.SQSTM1.5q = anno_simple(x = top_col_anno_df$CN.SQSTM1,
                                 #                            simple_anno_size = unit(3, "mm"), 
                                 #                            col = colors_cn),
                                 # Mut.VHL = anno_simple(x = top_col_anno_df$Mut.VHL,
                                 #                       simple_anno_size = unit(3, "mm"),
                                 #                       col = variant_class_colors),
                                 # Mut.PBRM1 = anno_simple(x = top_col_anno_df$Mut.PBRM1,
                                 #                         simple_anno_size = unit(3, "mm"),
                                 #                         col = variant_class_colors),
                                 # Mut.BAP1 = anno_simple(x = top_col_anno_df$Mut.BAP1,
                                 #                        simple_anno_size = unit(3, "mm"),
                                 #                        col = variant_class_colors),
                                 # Mut.PIK3CA = anno_simple(x = top_col_anno_df$Mut.PIK3CA,
                                 #                          simple_anno_size = unit(3, "mm"),
                                 #                          col = variant_class_colors),
                                 Cabozantinib = anno_simple(x = is_cab_treated,
                                                            simple_anno_size = unit(3, "mm"),
                                                            gp = gpar(col = "black"),
                                                            col = colors_is_cab_treated),
                                 Sapanisertib = anno_simple(x = is_sap_treated,
                                                            simple_anno_size = unit(3, "mm"),
                                                            gp = gpar(col = "black"),
                                                            col = colors_is_sap_treated),
                                 ACF = anno_simple(x = is_acf_treated,
                                                   simple_anno_size = unit(3, "mm"),
                                                   gp = gpar(col = "black"),
                                                   col = colors_is_acf_treated),
                                 Entinostat = anno_simple(x = is_ent_treated,
                                                          simple_anno_size = unit(3, "mm"),
                                                          gp = gpar(col = "black"),
                                                          col = colors_is_ent_treated),
                                 Control = anno_simple(x = is_control,
                                                       simple_anno_size = unit(3, "mm"),
                                                       gp = gpar(col = "black"),
                                                       col = colors_is_control),
                                 Bulk.WES = anno_simple(x = top_col_anno_df$WES,
                                                        simple_anno_size = unit(5, "mm"),
                                                        gp = gpar(col = "black"),
                                                        col = colors_seq_status),
                                 Bulk.RNA = anno_simple(x = top_col_anno_df$RNA,
                                                        simple_anno_size = unit(5, "mm"),
                                                        gp = gpar(col = "black"),
                                                        col = colors_seq_status),
                                 snRNA = anno_simple(x = top_col_anno_df$snRNA,
                                                     simple_anno_size = unit(5, "mm"),
                                                     gp = gpar(col = "black"),
                                                     col = colors_seq_status))

# make column split -------------------------------------------------------
col_split_vec <- model_ids
col_split_factor <- factor(x = col_split_vec, levels = c("RESL3", "RESL4", "RESL5", "RESL6", "RESL8", "RESL9", "RESL10", "RESL11", "RESL12"))
col_split_factor

# plot heatmap body with white-yellow-red ------------------------------------------------------
## make heatmap
p <- ComplexHeatmap::Heatmap(cluster_rows = F,
                             show_row_dend = F,
                             row_names_gp = grid::gpar(fontsize = row_fontsize),
                             column_split = col_split_factor, 
                             column_title_rot = 90, column_title_gp = grid::gpar(fontsize = 10),
                             top_annotation = top_col_anno,
                             show_column_names = F,
                             cluster_columns = F, 
                             show_column_dend = F,
                             matrix = plot_data_mat)
p
## make legend for top annotation
annotation_lgd = list(
  Legend(labels = names(colors_seq_status), 
         title = "Sequencing Status", 
         legend_gp = gpar(fill = colors_seq_status)))
## save heatmap
png(filename = paste0(dir_out, "test.", run_id, ".png"), 
    width = 1600, height = 800, res = 150)
### combine heatmap and heatmap legend
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = annotation_lgd)
dev.off()

