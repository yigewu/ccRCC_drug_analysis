# Yige Wu @WashU Apr 2020
## run clustering by for loop

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/load_pkgs.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/functions.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/variables.R")
library(Seurat)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
dir_out_elbowplot <- paste0(dir_out, "elbowplot", "/")
dir.create(dir_out_elbowplot)
dir_out_umap <- paste0(dir_out, "umap", "/")
dir.create(dir_out_umap)
# input dependencies ------------------------------------------------------
## input the path to the seurat-filtered srat objects
path_srat_obj_df <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/Resources/Analysis_Results/snrna_processing/sctransform/run_sctransform_by_for_loop/20200408.v1/path_seurat_sctransform_RDS.20200408.v1.tsv", data.table = F)

# process each sample -----------------------------------------------------
path_outputs <- NULL
for (sample_id in path_srat_obj_df$sample_id) {
  ## input filtered srat object
  path_srat <- path_srat_obj_df$path_output[path_srat_obj_df$sample_id == sample_id]
  srat <- readRDS(file = path_srat)
  
  # Run PCA
  srat <- RunPCA(object = srat, npcs = num_pcs)
  
  # Plot ans save the elbow plot, see if the number of PCs are reasonable
  p <- ElbowPlot(object = srat, ndims = num_pcs)
  file2write <- paste0(dir_out_elbowplot, sample_id, "_pc_elbowplot", run_id, ".png")
  png(filename = file2write, width = 800, height = 800, res = 150)
  print(p)
  dev.off()
  
  # Run UMAP
  srat <- RunUMAP(srat, dims = 1:num_pcs, reduction = "pca")
  
  # Determine the K-nearest neighbor graph
  srat <- FindNeighbors(object = srat, dims = 1:num_pcs)
  
  # Determine the clusters for various resolutions                                
  srat <- FindClusters(object = srat, resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
  
  for (res_tmp in c(0.4, 0.6, 0.8, 1.0, 1.4)) {
    # Assign identity of clusters
    Idents(object = srat) <- paste0("SCT_snn_res.", res_tmp)
    
    # Plot ans save the UMAP plot
    p <- DimPlot(srat,
                 reduction = "umap",
                 label = TRUE,
                 label.size = 6)
    
    file2write <- paste0(dir_out_umap, sample_id, "_", paste0("SCT_snn_res.", res_tmp), ".umap.", run_id, ".png")
    png(filename = file2write, width = 1000, height = 800, res = 150)
    print(p)
    dev.off()
  }
  
  ## save outputs
  file2write <- paste0(dir_out, sample_id, ".seruat_clustered.", run_id, ".RDS")
  saveRDS(object = srat, file = file2write, compress = T)
  
  ## store path to the outputs
  path_outputs <- c(path_outputs, file2write)
}
# make and write path to the quality metrics ------------------------------
path_outputs_df <- data.frame(sample_id = path_srat_obj_df$sample_id, path_output = path_outputs)
file2write <- paste0(dir_out, "path_seurat_clustered_RDS.", run_id, ".tsv")
write.table(x = path_outputs_df, file = file2write, quote = F, sep = "\t", row.names = F)

