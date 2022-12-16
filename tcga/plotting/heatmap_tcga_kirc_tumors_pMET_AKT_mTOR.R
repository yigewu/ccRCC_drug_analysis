# Yige Wu @ WashU 2019 Nov
## plot a heatmap with proteomics data from the discovery set data freeze


# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2",
  "ComplexHeatmap",
  "circlize"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
source("./ccRCC_drug_analysis//functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input data --------------------------------------------
expression_headers <- fread(data.table = F, 
                       input = "./Resources/Knowledge/gdac.broadinstitute.org_KIRC.Merge_protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data.Level_3.2016012800.0.0/KIRC.protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data.data.txt")

expression_df <- fread(data.table = F, skip = 1,
                       input = "./Resources/Knowledge/gdac.broadinstitute.org_KIRC.Merge_protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data.Level_3.2016012800.0.0/KIRC.protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data.data.txt")
colnames(expression_df) <- colnames(expression_headers)
rm(expression_headers)
# test_df <- fread(data.table = F, input = "~/Downloads/gdac.broadinstitute.org_KIRC.Merge_Clinical.Level_1.2016012800.0.0/KIRC.clin.merged.txt")
# test_df <- fread(data.table = F, input = "~/Downloads/gdac.broadinstitute.org_KIRC.Merge_Clinical.Level_1.2016012800.0.0/KIRC.merged_only_biospecimen_clin_format.txt")

# set parameters ----------------------------------------------------------
probes_plot <- c("c-Met_pY1235", "Akt_pS473", "Akt_pT308", "S6_pS235_S236", "S6_pS240_S244", "mTOR_pS2448", "p70S6K_pT389")
anno_fontsize <- 12

# make matrix for heatmap body --------------------------------------------
ncol(expression_df)
plot_data_df <- expression_df %>%
  filter(`Sample REF` %in% probes_plot)
plot_plot_mat <- as.matrix(plot_data_df[,-1])
rownames(plot_plot_mat) <- plot_data_df$`Sample REF`
## get column and row ids
column_ids_vec <- colnames(plot_plot_mat)

# make column annotation --------------------------------------------------
## calculate p_MET status
ca_data_df <- data.frame(sample_id = column_ids_vec,
                         pMET_score = as.vector(plot_plot_mat["c-Met_pY1235",]),
                         pMET_status = ifelse(as.vector(plot_plot_mat["c-Met_pY1235",]) > 0, "high", "low"),
                         AKTmTOR_score = colMeans(plot_plot_mat[c("Akt_pS473", "Akt_pT308", "S6_pS235_S236", "S6_pS240_S244", "mTOR_pS2448", "p70S6K_pT389"),]))
ca_data_df <- ca_data_df %>%
  mutate(AKTmTOR_status = ifelse(AKTmTOR_score > 0, "Activated", "non_Act")) %>%
  arrange(pMET_status, AKTmTOR_status, desc(pMET_score), AKTmTOR_score)
table(ca_data_df[, c("pMET_status", "AKTmTOR_status")])
nrow(ca_data_df)
## order columns
column_ids_ordered_vec <- ca_data_df$sample_id
plot_plot_mat <- plot_plot_mat[probes_plot,column_ids_ordered_vec]
## make colors
colors_pMET = c("high" = "orange", "low" = "grey90")
colors_AKTmTOR = c("Activated" = "red", "non_Act" = "lightblue")

## make annotation object
ca = HeatmapAnnotation(p_MET = anno_simple(x = ca_data_df$pMET_status, col = colors_pMET[ca_data_df$pMET_status]),
                       AKTmTOR = anno_simple(x = ca_data_df$AKTmTOR_status, col = colors_AKTmTOR[ca_data_df$AKTmTOR_status]))


## get color corresponding to values
col_protein = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# plot heatmap ------------------------------------------------------------
p <- Heatmap(matrix = plot_plot_mat, 
             name = "ProExpr", cluster_rows = F,
             top_annotation = ca, cluster_columns = F,
             show_column_names = F)
p
list_lgd = list(
  Legend(labels = names(colors_pMET), 
         title = "pMET", 
         title_gp = gpar(fontsize = anno_fontsize, fontface = "bold"),
         labels_gp = gpar(fontsize = anno_fontsize),
         grid_width = unit(5, "mm"),
         legend_gp = gpar(fill = colors_pMET), nrow = 2),
  Legend(labels = names(colors_AKTmTOR), 
         title = "AKTmTOR", 
         title_gp = gpar(fontsize = anno_fontsize, fontface = "bold"),
         labels_gp = gpar(fontsize = anno_fontsize),
         grid_width = unit(5, "mm"),
         legend_gp = gpar(fill = colors_AKTmTOR), nrow = 2))
  

  file2write <- paste0(dir_out, "TCGA", ".pdf")
  pdf(file2write, width = 10, height = 2.5, useDingbats = F)
  draw(object = p,
       annotation_legend_side = "right", annotation_legend_list = list_lgd)
  dev.off()


