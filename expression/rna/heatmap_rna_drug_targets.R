# Yige Wu @ WashU 2020 Feb
## plot gene expression for druggable targets across PDX lines
## two versions: one without scaling across lines, one with

# set up libraries and output directory -----------------------------------
## set up working directory and source functions and load libraries
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/ccRCC_Drug_shared.R")
## set run id 
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies --------------------------------------------
## input sample info to add in information about the model id, passage and TumorTissue info
sample_info_df <- readxl::read_excel(path = "./PDX-Pilot/DataFreeze/Sample_info/sampleInfo.v3/sampleInfo.washU_b1-b8.pdmr.other.passed.v3.20200130.from_jay.xlsx")
## input RNA expression
rna_exp_df <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/Resources/Analysis_Results/expression/rna/extract_rcc_rna_from_merged/20200216.v1/RCC_PDX.GeneExp.TPM.SampleID.tsv", data.table = F)

## input genes to plot
met_related_genes <- c("HGF", "MET", "AXL")
vegfr_related_genes <- c("FLT1", "KDR", "FLT3", "FLT4", "NRP1", "NRP2", 
                         "VEGFA", "VEGFB", "VEGFC", "VEGFD", "VEGFE")
cabo_related_genes <- c("KIT", "RET", "NTRK2", "TEK")

### FGF reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3042641/
fgfr_related_genes <- c("FGF1", "FGF2",
                        "FGF4", "FGF5", "FGF6",
                        "FGF3", "FGF7", "FGF10", "FGF22",
                        "FGF8", "FGF17", "FGF18",
                        "FGF9", "FGF16", "FGF320",
                        "FGF11", "FGF12", "FGF13", "FGF14",
                        "FGF19", "FGF21", "FGF23",
                        "FGFR1", "FGFR2", "FGFR3", "FGFR4")
pdgfr_related_genes <- c("PDGFA", "PDGFB", "PDGFC", "PDGFD",
                         "PDGFRA", "PDGFRB")

ras_related_genes <- c("KRAS", "HRAS", "NRAS", 
                       "BRAF", "RAF1", "MAP2K1", "MAPK1", "MAPK3",
                       "PIK3CA", "PIK3CD", "AKT1", "GSK3B")

mtor_related_genes <- c("MTOR", "AKTS1", "DEPTOR", "RPTOR", "MLST8", "MAPKAP1", "PRR5", "RICTOR",
                        "RPS6KB1", "EIF4E", "EIF4G", "EEF2",
                        "EIF4EBP", "EEF2K")
### for Acriflavine
#### get HIF targets
tf_tab <- fread(input = "./Ding_Lab/Projects_Current/TP53_shared_data/resources/PPI/TF_interactions.txt", data.table = F)
hif_targets <- tf_tab$target_genesymbol[tf_tab$source_genesymbol %in% c("HIF1A", "EPAS1")]
hif_related_genes <- c("HIF1A", "EPAS1", hif_targets)
### for Panobinostat and Entinostat
#### reference: https://pubchem.ncbi.nlm.nih.gov/compound/Panobinostat#section=Mechanism-of-Action&fullscreen=true
hdac_related_genes <- c("HDAC1", "HDAC2","HDAC3", "HDAC8",
                        "HDAC4", "HDAC5", "HDAC6", "HDAC7", "HDAC9", "HDAC10",
                        "HDAC8", "HDAC11")
### for Losartan
losartan_related_genes <- c("AGT", "AGTR1")

## sum up the genes
genes2plot <- c(met_related_genes, vegfr_related_genes, other_rtk_genes, fgfr_related_genes, pdgfr_related_genes, ras_related_genes, mtor_related_genes,
                hdac_related_genes, hif_related_genes, losartan_related_genes)
genes2plot <- unique(genes2plot)

## set plotting metrics
num_nonna <- 40
row_fontsize <- 9

# prepare matrix to plot --------------------------------------------------
## make the matrix to plot the heatmap
rna_exp_df %>% head()
exp_tab2plot <- rna_exp_df %>%
  filter(Gene %in% genes2plot)

exp_mat2plot <- exp_tab2plot %>%
  select(-Gene)
exp_mat2plot <- as.matrix(exp_mat2plot)
rownames(exp_mat2plot) <- exp_tab2plot$Gene

mat2plot <- log2(exp_mat2plot+1)
mat2plot <- mat2plot[intersect(genes2plot, rownames(mat2plot)),]
## make column annotation
# ca = HeatmapAnnotation(Sample_Type = ifelse(colnames(mat2plot) %in% tumor_bulk_aliquot_ids2plot, "Tumor", "Normal"),
#                        col = list(Sample_Type = c("Tumor" = "red", "Normal" = "green")))
# ca = ComplexHeatmap::HeatmapAnnotation(foo = anno_text(month.name, location = 0.5, just = "center",
#                                                        gp = gpar(fill = rep(2:4, each = 4), col = "white", border = "black"),
#                                                        width = max_text_width(month.name)*1.2))

## get color corresponding to values
# col_rna <- circlize::colorRamp2(c(min(mat2plot, na.rm=T), mean(mat2plot, na.rm=T), max(mat2plot, na.rm=T)),c("blue", "white", "red"))
col_rna <- circlize::colorRamp2(c(quantile(mat2plot, 0.1, na.rm=T), quantile(mat2plot, 0.5, na.rm=T), quantile(mat2plot, 0.9, na.rm=T)),c("blue", "white", "red"))
# col_rna <- circlize::colorRamp2(c(0,4),c("white", "red"))

# make heatmap ------------------------------------------------------------
## plot heatmap for MET related genes
mat_met <- mat2plot[(rownames(mat2plot) %in% met_related_genes),]
p_met <- ComplexHeatmap::Heatmap(mat_met,
                                 row_title = "MET related",
                                 row_title_gp = grid::gpar(fontsize = 12),
                                 col = col_rna,
                                 name = "log2(TPM+1)", 
                                 row_names_gp = grid::gpar(fontsize = row_fontsize),
                                 # top_annotation = ca,
                                 cluster_columns = F,
                                 cluster_rows = T)
png(filename = paste0(dir_out, "RCC_PDX.MET_related_genes.GeneExp.", run_id, ".png"), width = 800, height = 300 + 40*nrow(p_met), res = 150)
print(p_met)
dev.off()

## plot heatmap for VEGFR related genes
mat_vegfr <- mat2plot[(rownames(mat2plot) %in% vegfr_related_genes),]
p_vegfr <- ComplexHeatmap::Heatmap(mat_vegfr,
                                   row_title = "VEGFR related",
                                   row_title_gp = grid::gpar(fontsize = 12),
                                   col = col_rna,
                                   name = "log2(TPM+1)", 
                                   row_names_gp = grid::gpar(fontsize = row_fontsize),
                                   # top_annotation = ca,
                                   cluster_columns = F,
                                   cluster_rows = T)
p_vegfr
png(filename = paste0(dir_out, "RCC_PDX.VEGFR_related_genes.GeneExp.", run_id, ".png"), width = 800, height = 300 + 40*nrow(p_vegfr), res = 150)
print(p_vegfr)
dev.off()

## plot heatmap for cabo related
mat_cabo <- mat2plot[(rownames(mat2plot) %in% cabo_related_genes),]
p_cabo <- ComplexHeatmap::Heatmap(mat_cabo,
                                  row_title = "Cabozantinib related",
                                  row_title_gp = grid::gpar(fontsize = 12),
                                  col = col_rna,
                                  name = "log2(TPM+1)", 
                                  row_names_gp = grid::gpar(fontsize = 12),
                                  # top_annotation = ca,
                                  cluster_columns = F,
                                  cluster_rows = T)
png(filename = paste0(dir_out, "RCC_PDX.cabo_related_genes.GeneExp.", run_id, ".png"), width = 800, height = 300 + 40*nrow(p_cabo), res = 150)
print(p_cabo)
dev.off()

## plot heatmap for Losartan related
mat_losartan <- mat2plot[(rownames(mat2plot) %in% losartan_related_genes),]
p_losartan <- ComplexHeatmap::Heatmap(mat_losartan,
                                      row_title = "Losartan related",
                                      row_title_gp = grid::gpar(fontsize = 12),
                                      col = col_rna,
                                      name = "log2(TPM+1)", 
                                      row_names_gp = grid::gpar(fontsize = 12),
                                      # top_annotation = ca,
                                      cluster_columns = F,
                                      cluster_rows = T)
png(filename = paste0(dir_out, "RCC_PDX.Losartan_related_genes.GeneExp.", run_id, ".png"), width = 800, height = 300 + 40*nrow(p_losartan), res = 150)
print(p_losartan)
dev.off()

## plot heatmap for HIF related
mat_hif <- mat2plot[(rownames(mat2plot) %in% hif_related_genes),]
p_hif <- ComplexHeatmap::Heatmap(mat_hif,
                                 row_title = "HIF related",
                                 row_title_gp = grid::gpar(fontsize = 12),
                                 col = col_rna,
                                 name = "log2(TPM+1)", 
                                 row_names_gp = grid::gpar(fontsize = 12),
                                 # top_annotation = ca,
                                 cluster_columns = F,
                                 cluster_rows = T)
png(filename = paste0(dir_out, "RCC_PDX.HIF_related_genes.GeneExp.", run_id, ".png"), width = 800, height = 300 + 40*nrow(p_hif), res = 150)
print(p_hif)
dev.off()

## plot heatmap for HDAC
mat_hdac <- mat2plot[(rownames(mat2plot) %in% hdac_related_genes),]
p_hdac <- ComplexHeatmap::Heatmap(mat_hdac,
                                 row_title = "HDAC related",
                                 row_title_gp = grid::gpar(fontsize = 12),
                                 col = col_rna,
                                 name = "log2(TPM+1)", 
                                 row_names_gp = grid::gpar(fontsize = row_fontsize),
                                 # top_annotation = ca,
                                 cluster_columns = F,
                                 cluster_rows = T)
png(filename = paste0(dir_out, "RCC_PDX.HDAC_related_genes.GeneExp.", run_id, ".png"), width = 800, height = 300 + 40*nrow(p_hdac), res = 150)
print(p_hdac)
dev.off()

## plot heatmap for MTOR related genes
mat_mtor <- mat2plot[(rownames(mat2plot) %in% mtor_related_genes),]
p_mtor <- ComplexHeatmap::Heatmap(mat_mtor,
                                  row_title = "mTORC1/2 related",
                                  row_title_gp = grid::gpar(fontsize = 12),
                                  col = col_rna,
                                  name = "log2(TPM+1)", 
                                  row_names_gp = grid::gpar(fontsize = row_fontsize),
                                  # top_annotation = ca,
                                  cluster_columns = F,
                                  cluster_rows = T)
p_mtor
png(filename = paste0(dir_out, "RCC_PDX.MTOR_related_genes.GeneExp.", run_id, ".png"), width = 800, height = 300 + 40*nrow(p_mtor), res = 150)
print(p_mtor)
dev.off()

## plot heatmap for RAS related genes
mat_ras <- mat2plot[(rownames(mat2plot) %in% ras_related_genes),]
p_ras <- ComplexHeatmap::Heatmap(mat_ras,
                                   row_title = "RAS related",
                                   row_title_gp = grid::gpar(fontsize = 12),
                                   col = col_rna,
                                   name = "log2(TPM+1)", 
                                   row_names_gp = grid::gpar(fontsize = row_fontsize),
                                   # top_annotation = ca,
                                   cluster_columns = F,
                                   cluster_rows = T)
p_ras
png(filename = paste0(dir_out, "RCC_PDX.RAS_related_genes.GeneExp.", run_id, ".png"), width = 800, height = 60*nrow(mat_ras), res = 150)
print(p_ras)
dev.off()

## plot heatmap for PDGFR related genes
mat_pdgfr <- mat2plot[(rownames(mat2plot) %in% pdgfr_related_genes),]
p_pdgfr <- ComplexHeatmap::Heatmap(mat_pdgfr,
                                  row_title = "PDGFR related",
                                  row_title_gp = grid::gpar(fontsize = 12),
                                  col = col_rna,
                                  name = "log2(TPM+1)", 
                                  row_names_gp = grid::gpar(fontsize = row_fontsize),
                                  # top_annotation = ca,
                                  cluster_columns = F,
                                  cluster_rows = T)
p_pdgfr
png(filename = paste0(dir_out, "RCC_PDX.PDGFR_related_genes.GeneExp.", run_id, ".png"), width = 800, height = 80*nrow(p_pdgfr), res = 150)
print(p_pdgfr)
dev.off()



## plot heatmap for FGFR related genes
mat_fgfr <- mat2plot[(rownames(mat2plot) %in% fgfr_related_genes),]
p_fgfr <- ComplexHeatmap::Heatmap(mat_fgfr,
                                   row_title = "FGFR related",
                                   row_title_gp = grid::gpar(fontsize = 12),
                                   col = col_rna,
                                   name = "log2(TPM+1)", 
                                   row_names_gp = grid::gpar(fontsize = row_fontsize),
                                   # top_annotation = ca,
                                   cluster_columns = F,
                                   cluster_rows = T)
p_fgfr
png(filename = paste0(dir_out, "RCC_PDX.FGFR_related_genes.GeneExp.", run_id, ".png"), width = 800, height = 300 + 40*nrow(p_fgfr), res = 150)
print(p_fgfr)
dev.off()


dev.off()

