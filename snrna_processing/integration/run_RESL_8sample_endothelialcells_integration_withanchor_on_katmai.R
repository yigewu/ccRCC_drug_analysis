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
dir_rds <- "./Resources/Analysis_Results/snrna_processing/sctransform/run_sctransform_by_for_loop/20200416.v1/"
## input barcode2celltype annotation
barcode2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/map_barcode/map_barcode2celltype/20200501.v1/barcode2celltype_umapdata.20200501.v1.tsv")
id_model <- "RESL"
## set cell type to extract
celltype_short <- "Endothelial cells"
# get file names ----------------------------------------------------------
filenames_rds <- list.files(dir_rds) 
filenames_rds
filenames_rds <- filenames_rds[grepl(pattern = id_model, x = filenames_rds)]
filenames_rds

# input per object in for loop--------------------------------------------------------
list_srat <- list()
for (filename_rds in filenames_rds) {
  ## get sample id
  id_sample <- str_split_fixed(string = filename_rds, pattern = "\\.", 4)[,1] 
  id_sample
  
  ## get the barcode for the tumor cells
  barcode_tumorcells_df <- barcode2celltype_df %>%
    filter(Id_Sample == id_sample) %>%
    filter(Cell_Type.Short == celltype_short) %>%
    mutate(Barcode_Individual = str_split_fixed(string = Barcode_Integrated, pattern = '_', n = 2)[,1])
  
  ## get the path for RDS file
  path_rds <- paste0(dir_rds, filename_rds)
  
  ## read RDS file
  srat <- readRDS(file = path_rds)
  
  ## subset: tumor cells only
  srat_filtered <- subset(x = srat, cells = barcode_tumorcells_df$Barcode_Individual)
  
  ##s tore in the list
  list_srat[[id_sample]] <- srat_filtered
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
file2write <- paste0(dir_out, id_model, ".", gsub(x = celltype_short, pattern = "_", replacement = ""), ".", "integration.withanchor.", run_id, ".RDS")
saveRDS(object = srat_integrated, file = file2write, compress = T)
