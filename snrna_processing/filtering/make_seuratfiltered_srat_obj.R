# Yige Wu @WashU Apr 2020
## filter the Cell Ranger filtered barcode by the same threshold
## create the seurat object

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/load_pkgs.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/functions.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the paths to the quality metrics because it has multiplet info
path_qm_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/Resources/Analysis_Results/quality_control/filtering/make_quality_metrics/20200401.v1/path_quality_metrics_per_barcode.20200401.v1.tsv", data.table = F)
## set directory to the cell ranger outputs
dir_cellranger_outputs <- "./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/Resources/snRNA_Processed_Data/Cell_Ranger/outputs/GRCh38-3.0.0.premrna_and_mm10_premrna/"
## set sample ids
sample_ids <- list.files(path = dir_cellranger_outputs)
sample_ids
## set filtering threshold
mitoRatio_max <- 0.1
nCount_RNA_min <- 1000
nCount_RNA_max <- 80000
nFeature_RNA_min <- 200
nFeature_RNA_max <- 10000

# process each sample -----------------------------------------------------
path_outputs <- NULL
for (sample_id in sample_ids) {
  ## input quality metrics file
  path_qm <- path_qm_df$path_output[path_qm_df$sample_id == sample_id]
  qm_df <- fread(input = path_qm, data.table = F)
  
  ## filter by multiplet, mitoRatio, nCount_RNA, nFeature_RNA
  qm_filtered_df <- qm_df %>%
    filter(!is.na(barcode_raw))
  rownames(qm_filtered_df) <- qm_filtered_df$barcode
  
  ## Make matrix directory
  dir_matrix <- paste0(dir_cellranger_outputs, sample_id, "/outs/filtered_feature_bc_matrix/")
  
  ## Read in matrix
  input <- Seurat::Read10X(data.dir = dir_matrix)
  srat <- Seurat::CreateSeuratObject(counts = input,
                                     project = sample_id,
                                     min.features = nFeature_RNA_min,
                                     meta.data = qm_filtered_df)
  rm(input)
  
  filtered_srat <- subset(x = srat, 
                          subset= (!is.na(call)) & (call != "Multiplet") & 
                            (nCount_RNA >= nCount_RNA_min) & (nCount_RNA <= nCount_RNA_max) &
                            (nFeature_RNA >= nFeature_RNA_min) & (nFeature_RNA <= nFeature_RNA_max) &
                            (mitoRatio < mitoRatio_max))
  rm(srat)
  
  ## save outputs
  file2write <- paste0(dir_out, sample_id, ".seruat_filtered.", run_id, ".RDS")
  saveRDS(object = filtered_srat, file = file2write, compress = T)
  
  ## store path to the outputs
  path_outputs <- c(path_outputs, file2write)
}

# make and write path to the quality metrics ------------------------------
path_outputs_df <- data.frame(sample_id = sample_ids, path_output = path_outputs)
file2write <- paste0(dir_out, "path_seurat_filtered_RDS.", run_id, ".tsv")
write.table(x = path_outputs_df, file = file2write, quote = F, sep = "\t", row.names = F)
