# Yige Wu @ WashU 2020 Feb
## plot gene expression for druggable targets across PDX lines

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

# input dependencies --------------------------------------------
## input RNA expression
rna_exp_df <- fread("./Resources/Analysis_Results/expression/rna/generate_bulk_rna_table/20200612.v1/RCC_PDX.geneExp.20200612.v1.tsv", data.table = F)
## input te bulk genomics events
bulk_profile_df <- fread(input = "./Resources/Analysis_Results/bulk_data_summary/merge_selected_mutation_cnv/20200612.v1/RCC_PDX.SMG_Mutation.Known_Genes_CNV.20200612.v1.tsv", data.table = F)
## input all ccRCC samples data info
data_status_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_status/write_bulk_data_status/20200612.v1/RCC_PDX_Related_Samples.PDXNet_B1_9.HTAN_B1.Data_Status.20200612.v1.tsv")

# input genes to plot -----------------------------------------------------
ras_pathway_genes <- c("KRAS", "HRAS", "NRAS", 
                       "BRAF", "RAF1", "MAP2K1", "MAPK1", "MAPK3")
pi3k_pathway_genes <- c("PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PIK3R2",
                        "AKT1", "AKT2", "AKT3",
                        "GSK3B")
mtor_pathway_genes <- c("MTOR", "AKTS1", "DEPTOR", "RPTOR", "MLST8", "MAPKAP1", "PRR5", "RICTOR",
                        "RPS6KB1", "EIF4E", "EIF4G", "EEF2",
                        "EIF4EBP", "EEF2K")

## sum up the genes
genes2plot <- c(pi3k_pathway_genes, mtor_pathway_genes)
genes2plot <- unique(genes2plot)

## set plotting metrics
num_nonna <- 0
row_fontsize <- 9

# get ids in order --------------------------------------------------------
## sort data status
data_status_filtered_df <- data_status_df %>%
  filter(!is.na(Analysis_ID.RNA) & (Analysis_ID.RNA %in% colnames(rna_exp_df))) %>%
  arrange(ModelID, TreatmentLengthGroup, Treated.Cab, Treated.Sap, Treated.ACF, Treated.Entinostat)
analysis_ids <- data_status_filtered_df$Analysis_ID.RNA
analysis_ids
sample_ids <- mapvalues(x = analysis_ids, from = data_status_df$Analysis_ID.RNA, to = data_status_df$SampleID.AcrossDataType)
sample_ids
model_ids <- mapvalues(x = analysis_ids, from = data_status_df$Analysis_ID.RNA, to = data_status_df$ModelID)
model_ids
passages <- mapvalues(x = analysis_ids, from = data_status_df$Analysis_ID.RNA, to = data_status_df$NCI_Passage)
passages

# prepare matrix to plot --------------------------------------------------
## make the matrix to plot the heatmap
rna_exp_df %>% head()
exp_tab2plot <- rna_exp_df %>%
  filter(gene_symbol %in% genes2plot)
exp_mat2plot <- exp_tab2plot %>%
  select(-gene_symbol)
exp_mat2plot <- as.matrix(exp_mat2plot)
## add row names
rownames(exp_mat2plot) <- exp_tab2plot$gene_symbol
## sort by column
exp_mat2plot <- exp_mat2plot[,analysis_ids]
## transform the value
mat2plot <- log2(exp_mat2plot+1)
mat2plot <- mat2plot[intersect(genes2plot, rownames(mat2plot)),]

# make top column annotation --------------------------------------------------
## sort bulk profile table by sample ids
rownames(bulk_profile_df) <- bulk_profile_df$SampleID.AcrossDataType
bulk_profile_df <- bulk_profile_df[sample_ids,]
## get columns to use
colnames_colanno <- colnames(bulk_profile_df)[grepl(pattern = "Mut|CN|Treated", x = colnames(bulk_profile_df))]
colnames_colanno
colnames_colanno <- c(colnames_colanno, "TumorTissue", "TreatmentLengthGroup", "NCI_Passage")
## filter by columns
top_col_anno_df <- bulk_profile_df[,colnames_colanno]
### get unique color for each model
uniq_model_ids <- unique(model_ids)
uniq_model_colors <- Polychrome::dark.colors(n = 12)[1:length(uniq_model_ids)]
names(uniq_model_colors) <- uniq_model_ids
## get unique color for different passages
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
                                 # Mut.SETD2 = anno_simple(x = top_col_anno_df$Mut.SETD2,
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
                                 Control = anno_simple(x = is_control,
                                                       simple_anno_size = unit(3, "mm"),
                                                       col = colors_is_control))

# make row split ----------------------------------------------------------
rownames_mat <- rownames(mat2plot)
row_split_vec <- ifelse(rownames_mat %in% pi3k_pathway_genes, "PI3K Pathway", "MTOR Pathway")
row_split_factor <- factor(x = row_split_vec, levels = c("PI3K Pathway", "MTOR Pathway"))

# make heatmap ------------------------------------------------------------
col_rna <- circlize::colorRamp2(c(quantile(mat2plot, 0.1, na.rm=T), quantile(mat2plot, 0.5, na.rm=T), quantile(mat2plot, 0.9, na.rm=T)),c("blue", "white", "red"))
p <- ComplexHeatmap::Heatmap(mat2plot,
                             col = col_rna,
                             cluster_rows = T,
                             show_row_dend = F,
                             row_names_gp = grid::gpar(fontsize = row_fontsize),
                             row_split = row_split_factor,
                             row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 12),
                             column_split = col_split_factor, 
                             column_title_rot = 90, column_title_gp = grid::gpar(fontsize = 10),
                             top_annotation = top_col_anno,
                             show_column_names = F,
                             cluster_columns = F, 
                             show_column_dend = F,
                             name = "log2(TPM+1)")
p

png(filename = paste0(dir_out, "RCC_PDX.PI3K_MTOR_genes.GeneExp.", run_id, ".png"), width = 1500, height = 800, res = 150)
# print(p)
### combine heatmap and heatmap legend
draw(object = p)
dev.off()

