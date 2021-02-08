# Yige Wu @WashU Apr 2020
## for integrating 4 snRNA datasets for RESL5 on katmai

# set up libraries and output directory -----------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set time stamp for log file
timestamp <- paste0(run_id, ".", format(Sys.time(), "%H%M%S"))
## getting the path to the current script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
path_this_script <- thisFile()
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
library(Seurat)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)
## set log file
sink(file = paste0(dir_out, "Log.", timestamp, ".txt"))

# set dependencies --------------------------------------------------------
## set the directory containing the SCTransformed seurat objects in RDS file format
path_rds_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/filtering/make_mouse_seuratfiltered_srat_obj_katmai/20210208.v1/Path_to_Seurat_Objects.MouseCells.Filtered.20210208.v1.tsv")

# input per object in for loop--------------------------------------------------------
list_srat <- list()
for (id_sample_tmp in path_rds_df$id_sample) {
  ## get the path for RDS file
  path_rds <-path_rds_df$path_output_relative[path_rds_df$id_sample == id_sample_tmp]
  ## read RDS file and store in the list
  list_srat[[id_sample_tmp]] <- readRDS(file = path_rds)
  print(dim(list_srat[[id_sample_tmp]]))
}
cat("Finished creating the seurat object list!\n\n\n")
# start integration -------------------------------------------------------
# Select the most variable features to use for integration
features_integ <- SelectIntegrationFeatures(object.list = list_srat, 
                                            nfeatures = num_var_features, verbose = T) 
cat("Finished SelectIntegrationFeatures!\n\n\n")
# Prepare the SCT list object for integration
list_srat <- PrepSCTIntegration(object.list = list_srat, 
                                anchor.features = features_integ, verbose = T)
cat("Finished PrepSCTIntegration!\n\n\n")

# Find best buddies - can take a while to run
anchors_integ <- FindIntegrationAnchors(object.list = list_srat, 
                                        normalization.method = "SCT", 
                                        anchor.features = features_integ, verbose = T)
cat("Finished FindIntegrationAnchors!\n\n\n")
rm(list_srat)
# Integrate across conditions
srat <- IntegrateData(anchorset = anchors_integ, 
                      normalization.method = "SCT", verbose = T)
cat("Finished IntegrateData!\n\n\n")
rm(anchors_integ)

# dimension reduction -----------------------------------------------------
## PCA
### run PCA
srat <- RunPCA(object = srat)
cat("Finished RunPCA!\n\n\n")

## UMAP
### run UMAP
srat <- RunUMAP(srat, dims = 1:num_pcs, reduction = "pca")
### visualize UMAP
p <- DimPlot(object = srat, split.by = "orig.ident", reduction = "umap")
file2write <- paste0(dir_out, "umap.bysample", ".png")
png(filename = file2write, width = 4000, height = 1000, res = 150)
print(p)
dev.off()
cat("Finished RunUMAP!\n\n\n")

# cluster the cells -------------------------------------------------------
# Determine the K-nearest neighbor graph
srat <- FindNeighbors(object = srat, dims = 1:num_pcs)
cat("Finished FindNeighbors!\n\n\n")

# Determine the clusters for various resolutions                                
srat <- FindClusters(object = srat, resolution = findclusters_res)
cat("Finished FindClusters!\n\n\n")

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "MouseCells_8sample_integration.withanchor.", run_id, ".RDS")
saveRDS(object = srat, file = file2write, compress = T)
cat("Finished saveRDS!\n\n\n")
sink()
