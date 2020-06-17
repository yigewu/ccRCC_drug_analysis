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
## input te bulk genomics/methylation events
# bulk_profile_df <- fread(input = "./Resources/Analysis_Results/bulk_data_summary/merge_selected_mutation_cnv/20200526.v1/RCC_PDX.SMG_Mutation.Known_Genes_CNV.20200526.v1.tsv", data.table = F)
bulk_profile_df <- fread(input = "./Resources/Analysis_Results/bulk_data_summary/merge_selected_mutation_cnv/20200528.v1/RCC_PDX.SMG_Mutation.Known_Genes_CNV.20200528.v1.tsv", data.table = F)

# make data matrix for heatmap body ---------------------------------------
## filter samples into only drug treated ones
bulk_profile_df <- bulk_profile_df %>%
  filter(TumorTissue %in% c("PDX tumor;control", "PDX tumor;treated"))
## sort
bulk_profile_df <- bulk_profile_df %>%
  arrange(ModelID, NCI_Passage, Treated.Cab, Treated.Sap, Treated.ACF, Treated.Entinostat, TreatmentLengthGroup)
## get unique analysis ids
analysis_ids <- bulk_profile_df$Analysis_ID.WES
## make empty data frmae
plot_data_mat <- matrix(nc = length(analysis_ids), nr = 0)
plot_data_mat %>% head()

# make top column annotation --------------------------------------------------
## get columns to use
colnames_colanno <- colnames(bulk_profile_df)[grepl(pattern = "Mut|CN|Treated", x = colnames(bulk_profile_df))]
colnames_colanno
colnames_colanno <- c(colnames_colanno, "TumorTissue", "TreatmentLengthGroup", "NCI_Passage")
## filter by columns
top_col_anno_df <- bulk_profile_df[,colnames_colanno]
## add row names
rownames(top_col_anno_df) <- analysis_ids
### get model name
model_ids <- bulk_profile_df$ModelID
### get unique color for each model
uniq_model_ids <- unique(model_ids)
uniq_model_colors <- Polychrome::dark.colors(n = 12)[1:length(uniq_model_ids)]
names(uniq_model_colors) <- uniq_model_ids
## get unique color for different passages
passages <- top_col_anno_df$NCI_Passage
uniq_passages <- unique(passages)
uniq_passages
colors_uniq_passages <- RColorBrewer::brewer.pal(n = 5, name = "PuRd")[3:5]
names(colors_uniq_passages) <- c("P3", "P4", "P5")
## make color for VHL CNV
colors_cn <- rev(RColorBrewer::brewer.pal(n = 9, name = "PuOr"))[3:9]
names(colors_cn) <- 0:6
## make color for treatment length
treatmentlength <- top_col_anno_df$TreatmentLengthGroup
colors_treatmentlength <- RColorBrewer::brewer.pal(n = 4, name = "BuPu")[3:4]
names(colors_treatmentlength) <- c("2_Months", "1_Month")
## make color for treatment info
colors_treatment <- RColorBrewer::brewer.pal(n = 5, name = "Set2")
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
top_col_anno = HeatmapAnnotation(ModelID = anno_text(x = model_ids, 
                                                     location = 0.5, just = "center",
                                                     gp = gpar(fill = uniq_model_colors[model_ids], col = "white", border = "black"),
                                                     width = max_text_width(model_ids)*1.2),
                                 NCI_Passage = anno_text(x = passages, 
                                                         location = 0.5, just = "center",
                                                         gp = gpar(fill = colors_uniq_passages[passages], col = "white", border = "black"),
                                                         width = max_text_width(model_ids)*1.2),
                                 CN.VHL.3p = anno_simple(x = top_col_anno_df$CN.VHL,
                                                         simple_anno_size = unit(3, "mm"), 
                                                         col = colors_cn),
                                 CN.HIF1A.14q = anno_simple(x = top_col_anno_df$CN.HIF1A,
                                                            simple_anno_size = unit(3, "mm"), 
                                                            col = colors_cn),
                                 CN.SQSTM1.5q = anno_simple(x = top_col_anno_df$CN.SQSTM1,
                                                            simple_anno_size = unit(3, "mm"), 
                                                            col = colors_cn),
                                 Mut.VHL = anno_simple(x = top_col_anno_df$Mut.VHL,
                                                       simple_anno_size = unit(3, "mm"),
                                                       col = variant_class_colors),
                                 # Mut.PBRM1 = anno_simple(x = top_col_anno_df$Mut.PBRM1,
                                 #                         simple_anno_size = unit(3, "mm"),
                                 #                         col = variant_class_colors),
                                 # Mut.BAP1 = anno_simple(x = top_col_anno_df$Mut.BAP1,
                                 #                        simple_anno_size = unit(3, "mm"),
                                 #                        col = variant_class_colors),
                                 Mut.PIK3CA = anno_simple(x = top_col_anno_df$Mut.PIK3CA,
                                                          simple_anno_size = unit(3, "mm"),
                                                          col = variant_class_colors),
                                 TreatmentLength = anno_text(x = treatmentlength, 
                                                             location = 0.5, just = "center",
                                                             gp = gpar(fill = colors_treatmentlength[treatmentlength], col = "white", border = "black", fontsize = 8),
                                                             width = max_text_width(treatmentlength)*1.4),
                                 Cabozantinib = anno_simple(x = is_cab_treated,
                                                            simple_anno_size = unit(3, "mm"),
                                                            col = colors_is_cab_treated),
                                 Sapanisertib = anno_simple(x = is_sap_treated,
                                                            simple_anno_size = unit(3, "mm"),
                                                            col = colors_is_sap_treated),
                                 ACF = anno_simple(x = is_acf_treated,
                                                   simple_anno_size = unit(3, "mm"),
                                                   col = colors_is_acf_treated),
                                 Entinostat = anno_simple(x = is_ent_treated,
                                                          simple_anno_size = unit(3, "mm"),
                                                          col = colors_is_ent_treated),
                                 Is.Control = anno_simple(x = is_control,
                                                          simple_anno_size = unit(3, "mm"),
                                                          col = colors_is_control))


# plot heatmap body with white-yellow-red ------------------------------------------------------
## make heatmap
p <- Heatmap(matrix = plot_data_mat,
             top_annotation = top_col_anno,
             show_heatmap_legend = F)
p
## make legend for top annotation
annotation_lgd = list(
  Legend(labels = names(colors_cn), 
         title = "Bulk WES Copy Number", 
         legend_gp = gpar(fill = colors_cn)),
  Legend(labels = names(variant_class_colors), 
         title = "Bulk Mutation Class", 
         legend_gp = gpar(fill = variant_class_colors)))
## save heatmap
png(filename = paste0(dir_out, "test.", run_id, ".png"), 
    width = 1600, height = 800, res = 150)
### combine heatmap and heatmap legend
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = annotation_lgd)
dev.off()

