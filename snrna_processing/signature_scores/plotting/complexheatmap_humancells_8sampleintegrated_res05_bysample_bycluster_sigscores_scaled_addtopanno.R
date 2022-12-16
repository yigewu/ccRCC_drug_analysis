# Yige Wu @WashU May 2022

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "circlize",
  "ComplexHeatmap"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

# prepare for plotting ----------------------------------------------------
fontsize_plot <- 32
anno_fontsize <- 32
cellwidth <- 3
source("./ccRCC_drug_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input median(?) signature scores per cluster
results_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/signature_scores/other/make_median_signaturescores_bygeneset_bycluster_bysample/20220918.v1/median_scores.res05.bycluster.bysample.humancells.8sampleintegrated.20220918.v1.tsv")
## input the auto-correlation results
sigCorr_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/signature_scores/run_vision/run_vision_on_humancells_8sample_integrated/20220907.v1/humancells.8sampleintegrated.Vision.SignatureAutocorrelation.20220907.v1.tsv")
## input the annotation for the hallmark gene sets
hallmark_anno_df <- readxl::read_xlsx(path = "../ccRCC_snRNA/Resources/Knowledge/Databases/MSigDB/Hallmark_gene_sets_summary.xlsx")
##
sig_bycluster_long_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/signature_scores/compare_scores_across_clusters/wilcox_eachcluster_vs_rest_bygeneset_byres05/20220908.v1/wilcox.eachcluster_vs_rest.res1.humancells.8sampleintegrated.20220908.v1.tsv")

# make matrix data for heatmap body color-----------------------------------------------------------------
## get gene sets to plot
genesets_plot <- sigCorr_df$gene_set[sigCorr_df$FDR < 0.1 & sigCorr_df$C > 0.1]
sigCorr_filtered_df <- sigCorr_df %>%
  filter(FDR < 0.1 & C > 0.1) %>%
  arrange(desc(C))
# genesets_plot <- sigCorr_df$gene_set[sigCorr_df$FDR < 0.1 & sigCorr_df$C > 0.15]
genesets_plot <- genesets_plot[!grepl(pattern = "GOBP_", x = genesets_plot)]
# genesets_plot <- genesets_plot[!(genesets_plot %in% c("HALLMARK_UV_RESPONSE"))]
rm(sigCorr_df)
## extract the data for the matrix
plotdata_long_df <- results_df %>%
  filter(gene_set %in% genesets_plot) %>%
  mutate(cluster = str_split_fixed(string = group, pattern = "\\-", n = 4)[,4]) %>%
  mutate(model = str_split_fixed(string = group, pattern = "\\-", n = 4)[,1]) %>%
  mutate(treatment = str_split_fixed(string = group, pattern = "\\-", n = 4)[,3]) %>%
  mutate(treatment = gsub(x = treatment, pattern = '[0-9]', replacement = ""))
plotdata_long_df$model[plotdata_long_df$model == "RESL5E"] <- "RESL5"
plotdata_long_df$model[plotdata_long_df$model == "RESL10F"] <- "RESL10"
plotdata_long_df$treatment[plotdata_long_df$treatment == "Cabo"] <- "Cab"
plotdata_long_df$treatment[plotdata_long_df$treatment == "Cabo_Sap"] <- "CabSap"
plotdata_long_df <- plotdata_long_df %>%
  mutate(x_plot = gene_set) %>%
  mutate(cluster = as.numeric(cluster)) %>%
  arrange(cluster) %>%
  mutate(y_plot = paste0(model, "_", treatment, "_", paste0("MC", (cluster+1))))
unique(plotdata_long_df$y_plot)
rm(results_df)
plotdata_wide_df <- dcast(data = plotdata_long_df, formula = x_plot ~ y_plot, value.var = "value")
plotdata_raw_mat <- as.matrix(plotdata_wide_df[,-1])
## scale across clusters
plotdata_mat <- apply(plotdata_raw_mat, 1, scale)
rownames(plotdata_mat) <- colnames(plotdata_raw_mat)
colnames(plotdata_mat) <- plotdata_wide_df$x_plot
rm(plotdata_raw_mat)
# ## correlation between gene sets
# cor_bygeneset_mat <- cor(plotdata_mat, method = "spearman")
# Heatmap(matrix = cor_bygeneset_mat,
#         col = colors_correlation,
#         show_row_names = T, cluster_rows = T, show_row_dend = F,
#         ## column
#         show_column_dend = F,
#         cluster_columns = T,  show_column_names = T,
#         # column_names_gp = gpar(fontsize = fontsize_plot), 
#         # column_dend_side = "bottom",
#         # name = "Correlation",
#         show_heatmap_legend = F, use_raster = T)

## correlation between cell groups
plotdata_mat1 <- cor(t(plotdata_mat), method = "spearman")
## reorder clusters
clusters_ordered <- unique(plotdata_long_df$y_plot)
plotdata_mat <- plotdata_mat[clusters_ordered,]
plotdata_mat1 <- plotdata_mat1[clusters_ordered, clusters_ordered]
row_ids <- rownames(plotdata_mat)
## get gene sets
column_ids <- colnames(plotdata_mat)
row_labels_vec <- clusters_ordered
rm(plotdata_long_df)
rm(plotdata_wide_df)
model_ids_vec <- str_split_fixed(string = row_ids, pattern = "_", n = 3)[,1]
treatment_vec <- str_split_fixed(string = row_ids, pattern = "_", n = 3)[,2]
clusters_vec <- str_split_fixed(string = row_ids, pattern = "_", n = 3)[,3]

# make column labels ------------------------------------------------------
column_id_df <- data.frame(column_id = column_ids, column_label1 = str_split_fixed(string = column_ids, pattern = "_", n = 2)[,2])
column_id_df$column_label2 <- column_id_df$column_label1
column_id_df$column_label2[column_id_df$column_label1 == "ACTIVATION_OF_THE_MRNA_UPON_BINDING_OF_THE_CAP_BINDING_COMPLEX_AND_EIFS_AND_SUBSEQUENT_BINDING_TO_43S"] <- "TRANSLATION_INITIATION_COMPLEX"
column_id_df$column_label2[column_id_df$column_label1 == "SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE"] <- "SRP-DEPENDENT_PROTEIN_TARGETING_TO_ER_MEMBRANE"
column_id_df$column_label2[column_id_df$column_label1 == "RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY"] <- "RESPONSE_OF_GCN2_TO_AMINO_ACID_DEFICIENCY"
column_id_df$column_label2 <- gsub(x = column_id_df$column_label2, pattern = "_", replacement = " ")
col_labels_vec <- column_id_df$column_label2
  
# make colors -------------------------------------------------------------
## make colors for the heatmap body
summary(as.vector(plotdata_mat))
colors_sig_zscore = circlize::colorRamp2(breaks = c(-2, 0, 2), 
                                         colors = c("purple", "black", "yellow"))
## 
colors_correlation <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red")) 
## colors for model
colors_model <- RColorBrewer::brewer.pal(n = 3, name = "Dark2")[c(1,2)]
names(colors_model) <- c("RESL5", "RESL10")
## colors for treatment
color_red <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[1]
color_green <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[3]
color_yellow <- RColorBrewer::brewer.pal(n = 7, name = "Set2")[6]
color_grey <- "grey50"
colors_treatment <- c("CT" = color_grey, "Cab" = color_red, "Sap" = color_green, "CabSap" = color_yellow)
## colors for clusters
clusters_uniq_vec <- unique(clusters_vec)
colors_bycluster <- Polychrome::palette36.colors(n = (length(clusters_uniq_vec)+1))[-2]
names(colors_bycluster) <- clusters_uniq_vec

# make row annotation -----------------------------------------------------
row_ha <- rowAnnotation(Model = anno_simple(x = model_ids_vec, col = colors_model[model_ids_vec], width = unit(4*cellwidth, "mm")),
                        Treatment = anno_simple(x = treatment_vec, col = colors_treatment[treatment_vec], width = unit(4*cellwidth, "mm")),
                        Cluster = anno_simple(x = clusters_vec, col = colors_bycluster[clusters_vec], width = unit(4*cellwidth, "mm")), 
                        annotation_name_side = "top", annotation_name_gp = gpar(fontsize = anno_fontsize))


# make column annotation --------------------------------------------------
## process signature score by cluster
sig_bycluster_long_df <- sig_bycluster_long_df %>%
  filter(gene_set %in% genesets_plot) %>%
  mutate(value = median_sigScore) %>%
  mutate(cluster_name = paste0("MC", (cluster+1)))
sig_bycluster_wide_df <- dcast(data = sig_bycluster_long_df, formula = gene_set ~ cluster_name, value.var = "value")
sig_bycluster_raw_mat <- as.matrix(sig_bycluster_wide_df[,-1])
sig_bycluster_scaled_mat <- t(apply(sig_bycluster_raw_mat, 1, scale))
colnames(sig_bycluster_scaled_mat) <- colnames(sig_bycluster_raw_mat)
rownames(sig_bycluster_scaled_mat) <- sig_bycluster_wide_df$gene_set
sig_bycluster_scaled_mat <- sig_bycluster_scaled_mat[column_ids, ]
col_ha <- HeatmapAnnotation(foo = sig_bycluster_scaled_mat, 
                            geneset_name = anno_text(col_labels_vec, gp = gpar(fontsize = 28)),
                            col = list(foo = colors_sig_zscore), 
                            show_legend = F, 
                            annotation_name_gp = gpar(fontsize = anno_fontsize), 
                            annotation_name_side = "left", simple_anno_size = unit(4*cellwidth, "mm"))

# plot version 2 --------------------------------------------------------------------
ht_opt$ROW_ANNO_PADDING = unit(0, "cm")
ht_opt$COLUMN_ANNO_PADDING = unit(5, "cm")
# ht_opt$COLUMN_ANNO_PADDING = unit(0, "cm")

p <- Heatmap(matrix = plotdata_mat1,
             col = colors_correlation,
             width = ncol(plotdata_mat1)*unit(cellwidth, "mm"),
             height = nrow(plotdata_mat1)*unit(cellwidth, "mm"),
             show_row_names = F, cluster_rows = T, show_row_dend = F,
             ## column
             show_column_dend = F,
             cluster_columns = T, column_labels = row_labels_vec, show_column_names = F,
             # column_names_gp = gpar(fontsize = fontsize_plot), 
             # column_dend_side = "bottom",
             # name = "Correlation",
             show_heatmap_legend = F, use_raster = T)
p <- p + Heatmap(matrix = plotdata_mat, 
                 col = colors_sig_zscore,
                 cell_fun = function(j, i, x, y, w, h, fill) {
                   # if (plotdata_mat[i,j] >= quantile(plotdata_mat[,j], 0.75)+1.5*IQR(plotdata_mat[,j])) {
                   #   grid.text("*", x, y)
                   # }
                   # if (plotdata_mat[i,j] >= quantile(plotdata_mat[i,], 0.75)+1.5*IQR(plotdata_mat[i,])) {
                   #   grid.rect(x = x, y = y, width = w, height = h,
                   #             gp = gpar(col = "red", fill = NA))
                   # }
                   if (plotdata_mat[i,j] >= quantile(plotdata_mat[i,], 0.75)+1.5*IQR(plotdata_mat[i,]) | plotdata_mat[i,j] >= quantile(plotdata_mat[,j], 0.75)+1.5*IQR(plotdata_mat[,j])) {
                     grid.rect(x = x, y = y, width = w, height = h,
                               gp = gpar(col = "red", fill = NA))
                   }
                 },
                 width = ncol(plotdata_mat)*unit(3*cellwidth, "mm"), 
                 height = nrow(plotdata_mat)*unit(cellwidth, "mm"),
                 ## column
                 top_annotation = col_ha,
                 cluster_columns = T, show_column_dend = T, column_dend_side = "bottom", column_dend_height = unit(5*cellwidth, "mm"),
                 cluster_column_slices = F,  show_column_names = F,
                 # column_labels = col_labels_vec, column_names_gp = gpar(fontsize = 25), column_names_side = "top",
                 ## row
                 cluster_rows = F,
                 show_row_names = F,
                 left_annotation = row_ha,
                 # row_names_side = "left", 
                 row_labels = row_labels_vec,
                 show_heatmap_legend = F, use_raster = T)

list_lgd = list(
  Legend(labels = names(colors_bycluster), 
         title = "Meta-cluster", 
         title_gp = gpar(fontsize = anno_fontsize, fontface = "bold"),
         labels_gp = gpar(fontsize = anno_fontsize),
         grid_width = unit(1, "cm"),
         legend_gp = gpar(fill = colors_bycluster), nrow = 2),
  Legend(col_fun = colors_sig_zscore, 
         title = "Signature z-score",
         title_gp = gpar(fontsize = anno_fontsize, fontface = "bold"),
         labels_gp = gpar(fontsize = anno_fontsize),
         legend_width = unit(8, "cm"), 
         grid_height = unit(1, "cm"),
         direction = "horizontal"),
  Legend(col_fun = colors_correlation, 
         title = "Correlation",
         title_gp = gpar(fontsize = anno_fontsize, fontface = "bold"),
         labels_gp = gpar(fontsize = anno_fontsize),
         legend_width = unit(8, "cm"), 
         grid_height = unit(1, "cm"),
         direction = "horizontal"),
  Legend(labels = names(colors_model), 
         title = "Model", 
         legend_gp = gpar(fill = colors_model),
         title_gp = gpar(fontsize = anno_fontsize, fontface = "bold"),
         labels_gp = gpar(fontsize = anno_fontsize),
         grid_width = unit(1, "cm"),
         nrow = 1),
  Legend(labels = c("CT", "Cab", "Sap", "Cab+Sap"), 
         title = "Treatment", 
         legend_gp = gpar(fill = colors_treatment),
         title_gp = gpar(fontsize = anno_fontsize, fontface = "bold"),
         labels_gp = gpar(fontsize = anno_fontsize),
         grid_width = unit(1, "cm"),
         nrow = 1))


file2write <- paste0(dir_out, "heatmap2", ".pdf")
pdf(file2write, width = 40, height = 40, useDingbats = F)
draw(object = p,
     annotation_legend_side = "top", annotation_legend_list = list_lgd)
dev.off()
 
