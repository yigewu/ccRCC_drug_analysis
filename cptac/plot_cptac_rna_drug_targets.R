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

# input protein data ------------------------------------------------------
rna_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/mRNA/RNA_rpkm_tumor_normal.tsv", data.table = F)

# input bulk meta data ----------------------------------------------------
bulk_meta_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/cptac-metadata.csv")

# get the aliquot IDs for bulk corresponding to the snRNA aliquots --------
normal_bulk_aliquot_ids2plot <- bulk_meta_tab$RNA.ID[bulk_meta_tab$Set.A == "yes" & bulk_meta_tab$Type == "Normal" & !is.na(bulk_meta_tab$RNA.ID)]
normal_bulk_aliquot_ids2plot
case_ids2plot <- mapvalues(x = normal_bulk_aliquot_ids2plot, from = bulk_meta_tab$RNA.ID, to = bulk_meta_tab$Case.ID)
case_ids2plot
tumor_bulk_aliquot_ids2plot <- mapvalues(x = case_ids2plot, from = bulk_meta_tab$Case.ID[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"], to = as.vector(bulk_meta_tab$RNA.ID[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"]))
tumor_bulk_aliquot_ids2plot


# set NonNA threshold -----------------------------------------------------
num_nonna <- 40
row_fontsize <- 9

# input gene list ---------------------------------------------------------
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
genes2plot <- c(met_related_genes, vegfr_related_genes, cabo_related_genes, fgfr_related_genes, pdgfr_related_genes, ras_related_genes, mtor_related_genes,
                hdac_related_genes, hif_related_genes, losartan_related_genes)
genes2plot <- unique(genes2plot)

## make the matrix to plot the heatmap
exp_tab2plot <- rna_tab %>%
  filter(geneID %in% genes2plot) %>%
  select("geneID", tumor_bulk_aliquot_ids2plot, normal_bulk_aliquot_ids2plot)

exp_mat2plot <- exp_tab2plot %>%
  select(-geneID)
exp_mat2plot <- as.matrix(exp_mat2plot)
rownames(exp_mat2plot) <- exp_tab2plot$geneID

mat2plot <- log10(exp_mat2plot)
mat2plot <- mat2plot[intersect(genes2plot, rownames(mat2plot)),]
## make column annotation
ca = HeatmapAnnotation(Sample_Type = ifelse(colnames(mat2plot) %in% tumor_bulk_aliquot_ids2plot, "Tumor", "Normal"),
                       col = list(Sample_Type = c("Tumor" = "red", "Normal" = "green")))

## get color corresponding to values
# col_rna <- colorRamp2(c(min(mat2plot, na.rm=T), mean(mat2plot, na.rm=T), max(mat2plot, na.rm=T)),c("blue", "white", "red"))
# col_rna <- colorRamp2(c(quantile(mat2plot, 0.1, na.rm=T), mean(mat2plot, na.rm=T), quantile(mat2plot, 0.9, na.rm=T)),c("blue", "white", "red"))
# col_rna <- colorRamp2(c(-1.5, 0, 1.5),c("blue", "white", "red"))
col_rna <- colorRamp2(c(-2, 0, 2),c("blue", "white", "red"))

## plot heatmaps
p_hdac <- Heatmap(mat2plot[(rownames(mat2plot) %in% hdac_related_genes) & (rowSums(!is.na(mat2plot)) >= num_nonna),],
                      row_title = "HDACs",
                      row_title_gp = gpar(fontsize = 12),
                      col = col_rna,
                      name = "log10(RPKM)", 
                      row_names_gp = gpar(fontsize = row_fontsize),
                      top_annotation = ca,
                      cluster_columns = F,
                      cluster_rows = T)

p_hdac
## save heatmap to PNG
file2write <- paste0(dir_out, "HDAC_mRNA_Expression_log10RPKM.", run_id, ".png")
png(file2write, width = 2000, height = 600, res = 150)
print(p_hdac)
dev.off()


