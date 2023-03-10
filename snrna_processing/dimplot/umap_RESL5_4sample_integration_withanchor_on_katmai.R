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
## set the path to the rds file for integrated object
path_rds <- "./Resources/Analysis_Results/snrna_processing/clustering/cluster_RESL5_4sample_integration_withanchor_on_katmai/20200417.v1/RESL5_4sample_integration.withanchor.20200416.v1.clustered.RDS"
## input RDS file
srat <- readRDS(file = path_rds)

# make Dimplot, split ------------------------------------------------------------
p <- DimPlot(object = srat, reduction = "umap", label = T)
file2write <- paste0(dir_out, "umap.", "RESL5_4sample_integration.withanchor.", run_id, ".png")
png(filename = file2write, width = 2000, height = 1500, res = 150)
print(p)
dev.off()

# make Dimplot, split ------------------------------------------------------------
p <- DimPlot(object = srat, group.by = "call", split.by = "seurat_clusters", reduction = "umap", label = T)
file2write <- paste0(dir_out, "umap.", "RESL5_4sample_integration.withanchor.split.", run_id, ".png")
png(filename = file2write, width = 4000, height = 1000, res = 150)
print(p)
dev.off()


