# Yige Wu @WashU Apr 2020
## for making dimplot for RESL5_4sample_integration

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
## set integration id
id_integration <- "RESL5_4sample_integration.withanchor.20200409.v1"
## set the path to the rds file for integrated object
path_rds <- "./Resources/Analysis_Results/snrna_processing/integration/run_RESL5_4sample_integration_withanchor_on_katmai/20200409.v1/RESL5_4sample_integration.withanchor.20200409.v1.RDS"
## input RDS file
srat <- readRDS(file = path_rds)

# dimension reduction -----------------------------------------------------
## PCA
### run PCA
srat <- RunPCA(object = srat)
### visualize PCs
p <- PCAPlot(object = srat, group.by = "call", split.by = "orig.ident")
file2write <- paste0(dir_out, "pcaplot.", id_integration, ".", ".png")
png(filename = file2write, width = 4000, height = 1000, res = 150)
print(p)
dev.off()

## UMAP
### run UMAP
srat <- RunUMAP(srat, dims = 1:num_pcs, reduction = "pca")
### visualize UMAP
p <- DimPlot(object = srat, group.by = "call", split.by = "orig.ident", reduction = "umap")
file2write <- paste0(dir_out, "umap.", id_integration, ".", ".png")
png(filename = file2write, width = 4000, height = 1000, res = 150)
print(p)
dev.off()

# cluster the cells -------------------------------------------------------
# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:num_pcs)

# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))

# save output -------------------------------------------------------------


