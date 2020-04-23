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
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
source("./ccRCC_drug_analysis/plotting.R")
library(Seurat)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# set dependencies --------------------------------------------------------
## set integration id
id_integration <- "RESL5_4sample_integration.withanchor.20200417.v1"
## set the path to the rds file for integrated object
path_rds <- "./Resources/Analysis_Results/snrna_processing/clustering/cluster_RESL5_4sample_integration_withanchor_on_katmai/20200417.v1/RESL5_4sample_integration.withanchor.20200416.v1.clustered.RDS"
## input RDS file
srat <- readRDS(file = path_rds)

# fetch dats --------------------------------------------------------------
fetcheddata_df <- FetchData(object = srat, vars = c("call", "seurat_clusters", "UMAP_1", "UMAP_2"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, id_integration, ".fetched_data.", run_id, ".tsv")
write.table(x = fetcheddata_df, file = file2write, row.names = T, quote = F, sep = "\t")


