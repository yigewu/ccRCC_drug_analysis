# Yige Wu @WashU Apr 2020
## reference: https://satijalab.org/seurat/archive/v2.4/merge_vignette.html

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
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
# dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_Drug/"
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
path_rds_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/filtering/make_human_seuratfiltered_srat_obj_katmai/20210208.v1/Path_to_Seurat_Objects.HumanCellsss.Filtered.20210208.v1.tsv")

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

# start merging -------------------------------------------------------
## integrate without anchor
srat <- merge(x = list_srat[[1]], y = list_srat[2:length(list_srat)], project = "integrated")
rm(list_srat)
## scale data with all the features
srat <- SCTransform(srat, vars.to.regress = c("mitoRatio", 'nFeature_RNA', "nCount_RNA", 'S.Score', 'G2M.Score'), return.only.var.genes = F)
cat("Finished SCTransform!\n")
## keep it consistant with individual processing pipeline
srat <- RunPCA(srat, npcs = num_pcs, verbose = FALSE)
cat("Finished RUNPCA!\n")
srat <- RunUMAP(srat, reduction = "pca", dims = 1:num_pcs)
cat("Finished RUNUMAP!\n")
srat <- FindNeighbors(srat, reduction = "pca", dims = 1:num_pcs, force.recalc = T)
cat("Finished FindNeighbors!\n")
srat <- FindClusters(srat, resolution = findclusters_res)
cat("Finished FindClusters!\n")

# visualize-----------------------------------------------------
### visualize UMAP
p <- DimPlot(object = srat, split.by = "orig.ident", reduction = "umap")
file2write <- paste0(dir_out, "umap.bysample", ".png")
png(filename = file2write, width = 6000, height = 1000, res = 150)
print(p)
dev.off()
cat("Finished DimPlot!\n\n\n")

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "Humancells_8sample_merge.", run_id, ".RDS")
saveRDS(object = srat, file = file2write, compress = T)
cat("Finished saveRDS!\n\n\n")
sink()
