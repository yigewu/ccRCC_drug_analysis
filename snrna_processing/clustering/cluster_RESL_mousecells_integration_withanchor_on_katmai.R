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

# set dependencies --------------------------------------------------------
## set the path to the rds file for integrated object
path_rds <- "./Resources/Analysis_Results/snrna_processing/integration/run_RESL_8sample_mousecells_integration_withanchor_on_katmai/20200814.v1/RESL.mouse_cells.integration.withanchor.20200814.v1.RDS"
## set integration id
id_integration <- "RESL.mouse_cells.integration.withanchor.20200814.v1"
## input RDS file
srat <- readRDS(file = path_rds)

# dimension reduction -----------------------------------------------------
## PCA
### run PCA
srat <- RunPCA(object = srat)
### visualize PCs
p <- PCAPlot(object = srat, group.by = "call", split.by = "orig.ident", ncol = 4)
file2write <- paste0(dir_out, "pcaplot.", id_integration, ".png")
png(filename = file2write, width = 4000, height = 2000, res = 150)
print(p)
dev.off()

## UMAP
### run UMAP
srat <- RunUMAP(srat, dims = 1:num_pcs, reduction = "pca")
### visualize UMAP
p <- DimPlot(object = srat, group.by = "call", split.by = "orig.ident", reduction = "umap", ncol = 4)
file2write <- paste0(dir_out, "umap.", id_integration, ".png")
png(filename = file2write, width = 4000, height = 2000, res = 150)
print(p)
dev.off()

# cluster the cells -------------------------------------------------------
# Determine the K-nearest neighbor graph
srat <- FindNeighbors(object = srat, dims = 1:num_pcs)

# Determine the clusters for various resolutions                                
srat <- FindClusters(object = srat, resolution = findclusters_res)

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "clustered_srat.", id_integration, ".PC", num_pcs, ".Res", findclusters_res,  ".RDS")
saveRDS(object = srat, file = file2write, compress = T)

# fetch dats --------------------------------------------------------------
fetcheddata_df <- FetchData(object = srat, vars = c("orig.ident", "call", "seurat_clusters", "UMAP_1", "UMAP_2"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "umap_data.", id_integration, ".PC", num_pcs, ".Res", findclusters_res, ".tsv")
write.table(x = fetcheddata_df, file = file2write, row.names = T, quote = F, sep = "\t")


