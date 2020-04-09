# Yige Wu @ WashU 2019 Nov
## plot a heatmap with proteomics data from the discovery set data freeze

# set up libraries and output directory -----------------------------------
## set up working directory and source functions and load libraries
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/ccRCC_Drug_shared.R")
library(plyr)
library(ComplexHeatmap)
library(circlize)
## set run id 
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies --------------------------------------------
## input RNA
rna_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/mRNA/RNA_rpkm_tumor_normal.tsv", data.table = F)
## input protein
protein_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/proteome/6_CPTAC3_CCRCC_Whole_abundance_gene_protNorm=2_CB.tsv", data.table = F)
## input bulk meta data
bulk_meta_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/cptac-metadata.csv")
## set NonNA threshold
num_row_nonna <- 0
row_fontsize <- 9
## set genes to plot
genes2plot <- c("EZH1", "EZH2")

# plot RNA expression -----------------------------------------------------
## get the aliquot IDs for bulk RNA 
normal_bulk_rna_ids <- bulk_meta_tab$RNA.ID[bulk_meta_tab$Set.A == "yes" & bulk_meta_tab$Type == "Normal" & !is.na(bulk_meta_tab$RNA.ID)]
normal_bulk_rna_ids
case_ids2plot <- mapvalues(x = normal_bulk_rna_ids, from = bulk_meta_tab$RNA.ID, to = bulk_meta_tab$Case.ID)
case_ids2plot
tumor_bulk_rna_ids <- mapvalues(x = case_ids2plot, from = bulk_meta_tab$Case.ID[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"], to = as.vector(bulk_meta_tab$RNA.ID[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"]))
tumor_bulk_rna_ids
## make the matrix to plot the RNA expression
exp_tab2plot <- rna_tab %>%
  filter(geneID %in% genes2plot) %>%
  select("geneID", tumor_bulk_rna_ids, normal_bulk_rna_ids)

exp_mat2plot <- exp_tab2plot %>%
  select(-geneID)
exp_mat2plot <- as.matrix(exp_mat2plot)
rownames(exp_mat2plot) <- exp_tab2plot$geneID

mat2plot <- log10(exp_mat2plot+1)
mat2plot <- mat2plot[intersect(genes2plot, rownames(mat2plot)),]
## make column annotation
ca = HeatmapAnnotation(Sample_Type = ifelse(colnames(mat2plot) %in% tumor_bulk_rna_ids, "Tumor", "Normal"),
                       col = list(Sample_Type = c("Tumor" = "red", "Normal" = "green")))

## get color corresponding to values
col_rna <- colorRamp2(c(quantile(mat2plot, 0.1, na.rm=T), mean(mat2plot, na.rm=T), quantile(mat2plot, 0.9, na.rm=T)),c("blue", "white", "red"))

## plot heatmap for RNA expression
p <- Heatmap(mat2plot[(rowSums(!is.na(mat2plot)) >= num_row_nonna),],
             column_title = "EZH1/2 mRNA Expressin in CPTAC ccRCC Cohort",
             col = col_rna,
             name = "log10(RPKM+1)", 
             row_names_gp = gpar(fontsize = 9),
             top_annotation = ca,
             show_column_names = F,
             cluster_columns = F,
             cluster_rows = F)
## save plot for RNA expression
file2write <- paste0(dir_out, "EZH_mRNA_Expression_log10RPKM.", run_id, ".png")
png(file2write, width = 2000, height = 250, res = 150)
print(p)
dev.off()

# plot protein ------------------------------------------------------------
## get the aliquot IDs for bulk protein
normal_bulk_protein_ids <- bulk_meta_tab$Specimen.Label[bulk_meta_tab$Set.A == "yes" & bulk_meta_tab$Type == "Normal"]
normal_bulk_protein_ids
case_ids2plot <- mapvalues(x = normal_bulk_protein_ids, from = bulk_meta_tab$Specimen.Label, to = bulk_meta_tab$Case.ID)
case_ids2plot
tumor_bulk_protein_ids <- mapvalues(x = case_ids2plot, from = bulk_meta_tab$Case.ID[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"], to = as.vector(bulk_meta_tab$Specimen.Label[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"]))
tumor_bulk_protein_ids
## make the matrix to plot the heatmap
exp_tab2plot <- protein_tab %>%
  filter(Index %in% genes2plot) %>%
  select("Index", "ReferenceIntensity", tumor_bulk_protein_ids, normal_bulk_protein_ids)

protein_mat2plot <- exp_tab2plot %>%
  select(-Index) %>%
  select(-ReferenceIntensity)
rownames(protein_mat2plot) <- exp_tab2plot$Index

mat2plot <- as.matrix(protein_mat2plot) - as.vector(exp_tab2plot$ReferenceIntensity)
mat2plot <- mat2plot[intersect(genes2plot, rownames(mat2plot)),]
## make column annotation
ca = HeatmapAnnotation(Sample_Type = ifelse(colnames(mat2plot) %in% tumor_bulk_protein_ids, "Tumor", "Normal"),
                       col = list(Sample_Type = c("Tumor" = "red", "Normal" = "green")))

## get color corresponding to values
## get color corresponding to values
col_protein = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))

## plot heatmaps
p <- Heatmap(mat2plot[(rowSums(!is.na(mat2plot)) >= num_row_nonna),],
             column_title = "EZH1/2 Protein Expressin in CPTAC ccRCC Cohort",
             col = col_protein,
             name = "log2Intensity\n(Sample-Reference)", 
             row_names_gp = gpar(fontsize = 9),
             top_annotation = ca,
             show_column_names = F,
             cluster_columns = F,
             cluster_rows = F)
## save plot for RNA expression
file2write <- paste0(dir_out, "EZH_Protein_Expression_log10RPKM.", run_id, ".png")
png(file2write, width = 2000, height = 250, res = 150)
print(p)
dev.off()
