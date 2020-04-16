# Yige Wu @WashU Apr 2020
## filter the Cell Ranger filtered barcode by the same threshold
## create the seurat object

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the paths to the quality metrics because it has multiplet info
path_qm_df <- fread(input = "./Resources/Analysis_Results/snrna_processing/filtering/make_quality_metrics/20200401.v1/path_quality_metrics_per_barcode.20200401.v1.tsv", data.table = F)
### chanage the path
path_qm_df$path_output <- sapply(X = path_qm_df$path_output, FUN = function(path) {gsub(x = path, pattern = "quality_control", replacement = "snrna_processing")})
## set directory to the cell ranger outputs
dir_cellranger_outputs <- "./Resources/snRNA_Processed_Data/Cell_Ranger/outputs/GRCh38-3.0.0.premrna_and_mm10_premrna/"
## set sample ids
id_samples <- list.files(path = dir_cellranger_outputs)
id_samples
## set filtering threshold
mitoRatio_max <- 0.1
nCount_RNA_min_human <- 1000
nCount_RNA_min_mouse <- 500
nCount_RNA_max <- 80000
nFeature_RNA_min_human <- 200
nFeature_RNA_min_mouse <- 100
nFeature_RNA_max <- 10000

# process each sample -----------------------------------------------------
path_outputs <- NULL
summarisedqm_seuratfiltered_df <- NULL
for (id_sample in id_samples) {
  ## input quality metrics file
  path_qm <- path_qm_df$path_output[path_qm_df$sample_id == id_sample]
  qm_df <- fread(input = path_qm, data.table = F)
  
  ## filter by multiplet, mitoRatio, nCount_RNA, nFeature_RNA
  qm_cellrangerfiltered_df <- qm_df %>%
    filter(!is.na(barcode_raw))
  rownames(qm_cellrangerfiltered_df) <- qm_cellrangerfiltered_df$barcode
  ## Make matrix directory
  dir_matrix <- paste0(dir_cellranger_outputs, id_sample, "/outs/filtered_feature_bc_matrix/")
  
  ## Read in matrix
  input <- Seurat::Read10X(data.dir = dir_matrix)
  srat <- Seurat::CreateSeuratObject(counts = input,
                                     project = id_sample,
                                     meta.data = qm_cellrangerfiltered_df)
  rm(input)
  ## filter by quality metrics
  qm_seuratfiltered_df <- qm_cellrangerfiltered_df %>%
    filter(!is.na(call)) %>%
    filter(((call == "mm10_premrna") & (nCount_RNA >= nCount_RNA_min_mouse) & (nFeature_RNA >= nFeature_RNA_min_mouse)) | ((call == "GRCh38-3.0.0.premrna") & (nCount_RNA >= nCount_RNA_min_human) & (nFeature_RNA >= nFeature_RNA_min_human))) %>%
    filter((mitoRatio < mitoRatio_max) & (nCount_RNA <= nCount_RNA_max) & (nFeature_RNA <= nFeature_RNA_max))
  filtered_srat <- subset(x = srat, cells = qm_seuratfiltered_df$barcode)
  rm(srat)
  ## save outputs
  file2write <- paste0(dir_out, id_sample, ".seruat_filtered.", run_id, ".RDS")
  saveRDS(object = filtered_srat, file = file2write, compress = T)
  
  ## store path to the outputs
  path_outputs <- c(path_outputs, file2write)
  
  ## store the filtered quality metrics
  summarisedqm_seuratfiltered_tmp <- qm_seuratfiltered_df %>%
    group_by(call) %>%
    summarise(nBarcodes = n(), 
              nFeature_RNA_median = median(nFeature_RNA), 
              nCount_RNA_median = median(nCount_RNA)) %>%
    mutate(id_sample = id_sample)
  summarisedqm_seuratfiltered_df <- rbind(summarisedqm_seuratfiltered_df, summarisedqm_seuratfiltered_tmp)
}

# make and write path to the quality metrics ------------------------------
path_outputs_df <- data.frame(id_sample = id_samples, 
                              path_output_local = path_outputs, 
                              path_output_relative = gsub(x = path_outputs, pattern = dir_base, replacement = "./"))
file2write <- paste0(dir_out, "path_seurat_filtered_RDS.", run_id, ".tsv")
write.table(x = path_outputs_df, file = file2write, quote = F, sep = "\t", row.names = F)

# write the quality metrics after filtering ------------------------------
file2write <- paste0(dir_out, "summarised_quality_metrics.", run_id, ".tsv")
write.table(x = summarisedqm_seuratfiltered_df, file = file2write, quote = F, sep = "\t", row.names = F)



