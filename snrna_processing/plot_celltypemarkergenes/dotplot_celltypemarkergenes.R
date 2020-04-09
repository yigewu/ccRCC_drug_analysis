# Yige Wu @WashU Apr 2020
## make dotplot showing expression of cell type marker genes

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
source("./ccRCC_drug_analysis/plotting.R")
library(Seurat)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the path to the seurat-filtered srat objects
path_srat_obj_df <- fread("./Resources/Analysis_Results/snrna_processing/clustering/set_seurat_cluster_by_best_resolution/20200409.v1/path_seurat_clustered_best_reso_RDS.20200409.v1.tsv", data.table = F)
## input marker gene table
gene2celltype_df <- fread("./Resources/Gene_Lists/Cell_Type_Marker_Genes/RCC_Marker_Gene_Table - Gene2CellType_Tab.20200406.v1.tsv", data.table = F)

# make dotplot per sample -------------------------------------------------
for (id_sample in path_srat_obj_df$id_sample) {
  ## input srat object
  path_srat <- path_srat_obj_df$path_relative[path_srat_obj_df$id_sample == id_sample]
  srat <- readRDS(file = path_srat)
  
  ## set defailt assay to RNA
  DefaultAssay(srat) <- "RNA"
  
  ## get the genes within the cell type marker table
  genes2plot <- gene2celltype_df$Gene
  
  ## make feature names by adding prefix to the gene names
  features2plot <- paste0(prefix_human_feature, genes2plot)
  features2plot
  ## filter by feature names in the data
  features2plot <-  intersect(features2plot, srat@assays$RNA@counts@Dimnames[[1]])
  features2plot <- unique(features2plot)
  features2plot
  
  ## get the pct expressed for each gene in each cluster
  p <- DotPlot(object = srat, features = features2plot, col.min = 0)
  plot_data <- p$data
  ## transform the dataframe to matrix to better filter out genes with too low expressin
  plot_matrix <- dcast(data = plot_data, formula = features.plot ~ id, value.var = "pct.exp")
  plot_matrix %>% head()
  
  ## filter genes based on % expressed in each cluster
  ### for genes serving as markers for tumor cells, cutoff is 5%
  ### lower threshold than the rest because they are very important but they are lowly expressed
  genes2plot_tumor <- as.vector(gene2celltype_df$Gene[gene2celltype_df$Cell_Type_Group == "Malignant_Nephron_Epithelium"])
  features2plot_tumor <- paste0(prefix_human_feature, genes2plot_tumor)
  features_5pct <- plot_matrix[rowSums(plot_matrix[,unique(as.vector(plot_data$id))] > 5) >= 1, "features.plot"]
  features2plot_tumor_filtered <- intersect(features2plot_tumor, features_5pct)
  ### for other cell type marker genes, cutofff is 25%
  features_25pct <- plot_matrix[rowSums(plot_matrix[,unique(as.vector(plot_data$id))] > 25) >= 1, "features.plot"]
  features_25pct <- as.vector(features_25pct)
  ### merge filtered genes
  features2plot_filtered <- unique(c(features_25pct, features2plot_tumor_filtered))
  
  p <- DotPlot(object = srat, features = features2plot_filtered, col.min = 0)
  # p$data$gene_cell_type_group <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type_Group)
  # p$data$gene_cell_type1 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type1)
  # p$data$gene_cell_type2 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type2)
  # p$data$gene_cell_type3 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type3)
  # p$data$gene_cell_type4 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type4)
  p <- p  + RotatedAxis()
  # p <- p + facet_grid(.~gene_cell_type_group + gene_cell_type1 + gene_cell_type2 + gene_cell_type3 + gene_cell_type4, scales = "free", space = "free", drop = T)
  p <- p + theme(panel.spacing = unit(0, "lines"),
                 strip.background = element_blank(),
                 panel.border = element_rect(colour = "black"),
                 panel.grid.major = element_line(colour = "grey50"),
                 strip.text.x = element_text(angle = 0, vjust = 0.5),
                 axis.text.x = element_text(size = 15, face = "bold"),
                 strip.placement = "outside")
  p
  ## save plot
  file2write <- paste0(dir_out, id_sample, ".celltypemarkergene.dotplot,", run_id, ".png")
  png(file = file2write, width = 4000, height = 1200, res = 150)
  print(p)
  dev.off()
}

