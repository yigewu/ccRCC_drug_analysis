# Yige Wu @WashU Apr 2020
## set cluster number based on best resolution specific to each sample

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
library(Seurat)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the path to the seurat-filtered srat objects
path_srat_obj_df <- fread("./Resources/Analysis_Results/snrna_processing/clustering/run_clustering_by_for_loop/20200408.v1/path_seurat_clustered_RDS.20200408.v1.tsv", data.table = F)
## input the best resolution sheet
sample2res_df <- readxl::read_xlsx(path = "./Resources/snRNA_Processed_Data/Clustering/all_cells_clustering_different_resolution.xlsx")

# process each sample -----------------------------------------------------
path_outputs <- NULL
for (sample_id in path_srat_obj_df$sample_id) {
  ## input srat object
  path_srat <- path_srat_obj_df$path_output[path_srat_obj_df$sample_id == sample_id]
  srat <- readRDS(file = path_srat)
  stop("")
  ## get best resolution
  res_best <- sample2res_df$best_resolution[sample2res_df$sample_id == sample_id]
  res_best
  ## add best resolution cluster to meta data
  View(srat@meta.data)
  srat@meta.data$seurat_clusters <- srat@meta.data[, paste0("SCT_snn_res.", res_best)]
  Idents(srat) <- "seurat_clusters"

  ## save outputs
  file2write <- paste0(dir_out, sample_id, ".seurat_clustered_best_res.", run_id, ".RDS")
  saveRDS(object = srat, file = file2write, compress = T)
  
  ## store path to the outputs
  path_outputs <- c(path_outputs, file2write)
}
# make and write output path ------------------------------
path_outputs_df <- data.frame(id_sample = path_srat_obj_df$sample_id, 
                              path_full = path_outputs,
                              path_relative = gsub(x = path_outputs, pattern = dir_base, replacement = "./"))
file2write <- paste0(dir_out, "path_seurat_clustered_best_reso_RDS.", run_id, ".tsv")
write.table(x = path_outputs_df, file = file2write, quote = F, sep = "\t", row.names = F)

