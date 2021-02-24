# Yige Wu @WashU Feb 2021

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
library(Seurat)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# set dependencies --------------------------------------------------------
## set the directory containing the SCTransformed seurat objects in RDS file format
path_rds_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/filtering/make_human_seuratfiltered_srat_obj_katmai/20210208.v1/Path_to_Seurat_Objects.HumanCellsss.Filtered.20210208.v1.tsv")

# input per object in for loop--------------------------------------------------------
fetcheddata_sup_df <- NULL
for (id_sample_tmp in path_rds_df$id_sample) {
  ## get the path for RDS file
  path_rds <-path_rds_df$path_output_relative[path_rds_df$id_sample == id_sample_tmp]
  ## read RDS file
  srat <- readRDS(file = path_rds)
  ## fetch
  fetcheddata_df <- FetchData(object = srat, vars = c("orig.ident", "seurat_clusters", "UMAP_1", "UMAP_2"))
  fetcheddata_df$barcode <- rownames(fetcheddata_df)
  fetcheddata_sup_df <- rbind(fetcheddata_sup_df, fetcheddata_df)
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "HumanCells.Individual_sample.", "umap_data.", run_id, ".tsv")
write.table(x = fetcheddata_sup_df, file = file2write, row.names = T, quote = F, sep = "\t")


