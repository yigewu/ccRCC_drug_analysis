# Yige Wu @WashU Apr 2020
## for integrating 8 snRNA datasets on local
## making a template script transferable for running on cluster

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
## set the directory containing the SCTransformed seurat objects in RDS file format
dir_rds <- "./Resources/Analysis_Results/snrna_processing/sctransform/run_sctransform_by_for_loop/20200408.v1/"

# get file names ----------------------------------------------------------
filenames_rds <- list.files(dir_rds) 
filenames_rds
filenames_rds <- filenames_rds[grepl(pattern = "RESL10", x = filenames_rds)]
filenames_rds

# input per object in for loop--------------------------------------------------------
list_srat <- list()
for (filename_rds in filenames_rds) {
  ## get sample id
  id_sample <- str_split_fixed(string = filename_rds, pattern = "\\.", 4)[,1] 
  id_sample
  ## get the path for RDS file
  path_rds <- paste0(dir_rds, filename_rds)
  ## read RDS file and store in the list
  list_srat[[id_sample]] <- readRDS(file = path_rds)
}

# start integration -------------------------------------------------------
# Select the most variable features to use for integration
features_integ <- SelectIntegrationFeatures(object.list = list_srat, 
                                            nfeatures = num_var_features) 
# Prepare the SCT list object for integration
list_srat <- PrepSCTIntegration(object.list = list_srat, 
                                anchor.features = features_integ)
# Find best buddies - can take a while to run
anchors_integ <- FindIntegrationAnchors(object.list = list_srat, 
                                        normalization.method = "SCT", 
                                        anchor.features = features_integ, verbose = T)
rm(list_srat)
# Integrate across conditions
srat_integrated <- IntegrateData(anchorset = anchors_integ, 
                                   normalization.method = "SCT")

# Save integrated seurat object
file2write <- paste0(dir_out, "RESL10_4sample_integration.withanchor.", run_id, ".RDS")
saveRDS(object = srat_integrated, file = file2write, compress = T)
sink()