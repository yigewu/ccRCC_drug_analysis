# Yige Wu @WashU Apr 2020
## reference: https://satijalab.org/seurat/archive/v3.0/integration.html
## reference: https://hbctraining.github.io/scRNA-seq/lessons/06_SC_SCT_and_integration.html

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
path_rds_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/filtering/make_human_seuratfiltered_srat_obj_katmai/20210208.v1/Path_to_Seurat_Objects.HumanCellsss.Filtered.20210208.v1.tsv")

# input per object in for loop--------------------------------------------------------
list_srat <- list()
for (id_sample_tmp in path_rds_df$id_sample) {
  ## get the path for RDS file
  path_rds <-path_rds_df$path_output_relative[path_rds_df$id_sample == id_sample_tmp]
  ## read RDS file and store in the list
  srat <- readRDS(file = path_rds)
  cat(paste0("Finished readRDS for ", id_sample_tmp, " !\n\n\n"))
  print(dim(srat))
  
  ## get raw read count matrix
  raw_exp_mat <- seurat_object@assays$RNA@counts
  print(dim(raw_exp_mat))
  print(raw_exp_mat[1:5,1:4])

  ## write table
  file2write <- paste0(dir_out, id_sample_tmp, ".filtered.raw_counts.txt.gz")
  fwrite(x = raw_exp_mat, file = file2write, quote = F, sep = "\t", row.names = T)
  cat(paste0("Finished fwrite for ", id_sample_tmp, " !\n\n\n"))
}

sink()