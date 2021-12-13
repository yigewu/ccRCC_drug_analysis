# Yige Wu @ WashU 2021 Nov

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
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/0.sample_info/v3.20210116/RCC_PDX_Samples.20210115.v2.xlsx")
## input mutation
mutation_df <- fread(data.table = F, input = "./Resources/Analysis_Results/mutation/generate_bulk_mutation_table_for_data_freeze/20211122.v1/RCC_PDX.AAChange_VAF.20211122.v1.tsv")
## input CNV
cnv_arm_df <- fread(data.table = F, input = "./Data_Freeze/v1.dataFreeze.washU_rcc/2.copyNumber/gistic2.broad_values_by_arm.20200818.txt")
## input the protein meta data
protein_meta_data_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")

# process samples to be shown ---------------------------------------------
samples_plot_df <- sampleinfo_df %>%
  filter(ModelID %in% paste0("RESL", c(3, 4, 5, 10, 11, 12))) %>%
  # filter(DataType == "WES") %>%
  filter(Group == "PDX") %>%
  filter(ShortTag %in% c("Baseline", "Control", "Treated.Cabo", "Treated.Sap", "Treated.Cabo+Sap")) %>%
  mutate(column_id = paste0(ModelID, "_", ShortTag, "_", Treatment.Month))
columns_plot_df <- samples_plot_df %>%
  select(ModelID, ShortTag, Treatment.Month, column_id) %>%
  unique() %>%
  arrange(ModelID, Treatment.Month, factor(ShortTag, levels = c("Baseline", "Control", "Treated.Cabo", "Treated.Sap", "Treated.Cabo+Sap")))
protein_meta_data_df <- protein_meta_data_df %>%
  mutate(ModelID = str_split_fixed(string = `Sample ID`, pattern = "_", n = 3)[,1]) %>%
  mutate(Treatment.Month = str_split_fixed(string = `Treatment_length`, pattern = " ", n = 2)[,1]) %>%
  mutate(ShortTag = ifelse(Treatment %in% c("Cabo", "Sap"), paste0("Treated.", Treatment),
                           ifelse(Treatment == "Con", "Control", "Treated.Cabo+Sap"))) %>%
  mutate(column_id = paste0(ModelID, "_", ShortTag, "_", Treatment.Month)) %>%
  arrange(ModelID)
## get row/column ids
column_ids <- columns_plot_df$column_id
column_analysisids <- samples_plot_df$Analysis_ID

# make data matrix for heatmap body ---------------------------------------
## reformat data frame to matrix
plot_data_mat <- matrix(nc = length(column_ids), nr = 0)

# make colors -------------------------------------------------------------
## top column annotation object
color_gridline = "white"
color_na <- "grey50"
## make color for treatment type
RColorBrewer::display.brewer.pal(name = "Set1", n = 8)
colors_treatment <- c("#B2DF8A", RColorBrewer::brewer.pal(n = 5, name = "Set1")[c(1,2,4,3)])
names(colors_treatment) <- c("Baseline", "Treated.Cabo", "Treated.Sap", "Treated.Cabo+Sap", "Control")
## make color for treatment length
colors_treatmentlength <- c(RColorBrewer::brewer.pal(name = "Set1", n = 8)[c(7)], "#E6AB02", "#B2DF8A")
names(colors_treatmentlength) <- c("2 month", "1 month", "Baseline")
colors_mut = c("TRUE" = "#e7298a", "FALSE" = "grey90")
colors_cnstatus = c("Gain" = "#E41A1C", "Loss" = "#377EB8", "Neutral" = "white", "NA" = "grey90")
colors_datatype <- c("TRUE" = "#7570B3", "FALSE" = "grey90")

# make column annotation --------------------------------------------------
## make data status
wes_vec <- as.character(column_ids %in% samples_plot_df$column_id[samples_plot_df$DataType == "WES"])
rna_vec <- as.character(column_ids %in% samples_plot_df$column_id[samples_plot_df$DataType == "RNA-Seq"])
protein_vec <- as.character(column_ids %in% protein_meta_data_df$column_id)
snRNA_vec <- as.character(grepl(x = column_ids, pattern = "RESL5|RESL10") & grepl(x = column_ids, pattern = "_2"))
## make treatment status
treatment_vec <- mapvalues(x = column_ids, from = columns_plot_df$column_id, to = as.vector(columns_plot_df$ShortTag))
treatmentlength_vec <- mapvalues(x = column_ids, from = columns_plot_df$column_id, to = as.vector(columns_plot_df$Treatment.Month))
treatmentlength_vec <- paste0(treatmentlength_vec, " month"); treatmentlength_vec[treatmentlength_vec == "0 month"] <- "Baseline"

## make annotation object
top_col_anno = HeatmapAnnotation(
  Treatment = anno_simple(x = treatment_vec, simple_anno_size = unit(3, "mm"), col = colors_treatment),
  TreatmentLength = anno_simple(x = treatmentlength_vec, simple_anno_size = unit(3, "mm"), col = colors_treatmentlength[treatmentlength_vec]),
  gap1 = anno_empty(border = F, height = unit(0.5, "mm")),
  WES = anno_simple(x = wes_vec, simple_anno_size = unit(3.5, "mm"), col = colors_datatype),
  RNAseq = anno_simple(x = rna_vec, simple_anno_size = unit(3.5, "mm"), col = colors_datatype),
  Proteomics = anno_simple(x = protein_vec, simple_anno_size = unit(3.5, "mm"), col = colors_datatype),
  snRNASeq = anno_simple(x = snRNA_vec, simple_anno_size = unit(3.5, "mm"), col = colors_datatype),
  # gap2 = anno_empty(border = F, height = unit(0.5, "mm")),
  # Chr.3p = anno_simple(x = cnstatus_3p_vec, simple_anno_size = unit(3, "mm"), gp = gpar(col = color_gridline), col = colors_cnstatus),
  # Chr.5q = anno_simple(x = cnstatus_5q_vec, simple_anno_size = unit(3, "mm"), gp = gpar(col = color_gridline), col = colors_cnstatus),
  # Chr.14q = anno_simple(x = cnstatus_14q_vec, simple_anno_size = unit(3, "mm"), gp = gpar(col = color_gridline), col = colors_cnstatus),
  # gap3 = anno_empty(border = F, height = unit(0.5, "mm")),
  # VHL = anno_simple(x = vhl_mut_vec, simple_anno_size = unit(3, "mm"), gp = gpar(col = color_gridline), col = colors_mut),
  # PBRM1 = anno_simple(x = PBRM1_mut_vec, simple_anno_size = unit(3, "mm"), gp = gpar(col = color_gridline), col = colors_mut),
  # PIK3CA = anno_simple(x = PIK3CA_mut_vec, simple_anno_size = unit(3, "mm"), gp = gpar(col = color_gridline), col = colors_mut),
  # KDM5C = anno_simple(x = KDM5C_mut_vec, simple_anno_size = unit(3, "mm"), gp = gpar(col = color_gridline), col = colors_mut),
  # ARID1A = anno_simple(x = ARID1A_mut_vec, simple_anno_size = unit(3, "mm"), gp = gpar(col = color_gridline), col = colors_mut),
  annotation_name_side = "left")

# make column split -------------------------------------------------------
col_split_vec <- mapvalues(x = column_ids, from = columns_plot_df$column_id, to = as.vector(columns_plot_df$ModelID))
col_split_factor <- factor(x = col_split_vec, levels = paste0("RESL", c(3, 4, 5, 10, 11, 12)))

# make heatmap --------------------------------------------------------------------
p <- Heatmap(matrix = plot_data_mat,
             na_col = color_na,
             ## column
             column_split = col_split_factor, cluster_columns = F,
             column_title_gp = gpar(fontsize = 9), column_title_rot = 0,
             # bottom_annotation = bottom_col_anno, 
             top_annotation = top_col_anno,
             # column_gap = unit(x = 0, units = "mm"),
             ## row
             show_heatmap_legend = F)
p


# make legend -------------------------------------------------------------
annotation_lgd = list(
  Legend(labels = names(colors_treatment), 
         title = "Treatment", 
         legend_gp = gpar(fill = colors_treatment)),
  Legend(labels = names(colors_treatmentlength), 
         title = "Treatment\nlength", 
         legend_gp = gpar(fill = colors_treatmentlength)),
  Legend(labels = names(colors_datatype), 
         title = "Data\navailable", 
         legend_gp = gpar(fill = colors_datatype)))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "heatmap.png")
png(file2write, width = 1000, height = 600, res = 150)
draw(p, annotation_legend_side = "bottom", annotation_legend_list = annotation_lgd)  #Show the heatmap
dev.off()
file2write <- paste0(dir_out, "heatmap.pdf")
pdf(file2write, width = 6, height = 2.5, useDingbats = F)
draw(p, annotation_legend_side = "bottom", annotation_legend_list = annotation_lgd)  #Show the heatmap
dev.off()

