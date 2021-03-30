# Yige Wu @WashU Jan 2021
## filter the Cell Ranger filtered barcode by the same threshold
## create the seurat object

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
dir_out_elbowplot <- paste0(dir_out, "elbowplot", "/")
dir.create(dir_out_elbowplot)
dir_out_umap <- paste0(dir_out, "umap", "/")
dir.create(dir_out_umap)
## set log file
sink(file = paste0(dir_out, "Log.", timestamp, ".txt"))

# input dependencies ------------------------------------------------------
## set directory to the cell ranger outputs
dir_cellranger_outputs <- "./Resources/snRNA_Processed_Data/Cell_Ranger/cellranger-5.0.1_Ref-2020-A/"
## set sample ids
ids_sample <- list.files(path = dir_cellranger_outputs)
ids_sample <- ids_sample[grepl(pattern = 'RESL', x = ids_sample)]
# ids_sample
## input paths to the scrublet outputs
scrublet_out_paths_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/other/write_paths_to_scrublet_output_tables/20210127.v1/Paths_to_Scrulbet_Output_Tables.20210127.v1.tsv")
## input scrublet 
## set filtering threshold
mitoRatio_max <- 0.1
nCount_RNA_min_human <- 1000
nCount_RNA_min_mouse <- 500
nCount_RNA_max <- 80000
nFeature_RNA_min_human <- 200
nFeature_RNA_min_mouse <- 100
nFeature_RNA_max <- 10000

# process each sample -----------------------------------------------------
metadata_cellrangerfiltered_df <- NULL
path_outputs <- NULL
metadata_seuratfiltered_df <- NULL
for (id_sample in ids_sample) {
  print(id_sample)
  
  ## Make matrix directory
  dir_matrix <- paste0(dir_cellranger_outputs, id_sample, "/", id_sample, "/outs/filtered_feature_bc_matrix/")
  
  ## Read in matrix
  input <- Seurat::Read10X(data.dir = dir_matrix)
  srat <- Seurat::CreateSeuratObject(counts = input,
                                     project = id_sample)
  rm(input)
  
  ## Compute percent mito ratio
  srat$mitoRatio <- Seurat::PercentageFeatureSet(object = srat, pattern = "MT-")
  srat$mitoRatio <- srat@meta.data$mitoRatio / 100
  
  ## extract meta data
  metadata_df <- srat@meta.data
  metadata_df$barcode_raw <- rownames(metadata_df)
  metadata_df$barcode <- paste0(rownames(metadata_df), "-1")
  print(head(metadata_df))
  
  ## add barcode classification
  class_df <- fread(data.table = F, input = paste0(dir_cellranger_outputs, id_sample, "/", id_sample, "/outs/analysis/gem_classification.csv"))
  print(head(class_df))
  
  metadata_df <- merge(x = metadata_df, y = class_df, by.x = c("barcode"), by.y = c("barcode"), all.x = T)
  print(head(metadata_df))
  
  ## add scrublet information if any
  ### input scrublet output
  if (id_sample %in% scrublet_out_paths_df$Sample_id) {
    path_scrublet <- scrublet_out_paths_df$Path_relative[scrublet_out_paths_df$Sample_id == id_sample]
    scrublet_out_df <- fread(data.table = F, input = path_scrublet)
    metadata_df <- merge(x = metadata_df, y = scrublet_out_df, by.x = c("barcode"), by.y = c("Barcode"), all.x = T)
  } else {
    metadata_df <- metadata_df %>%
      mutate(doublet_score = NA) %>%
      mutate(predicted_doublet = FALSE)
  }
  print(head(metadata_df))
  ## merge into the super table
  metadata_cellrangerfiltered_df <- rbind(metadata_cellrangerfiltered_df, metadata_df)
  
  ## filter by quality metrics
  metadata_filtered_df1 <- metadata_df %>%
    filter((call == "mm10") & (nCount_RNA >= nCount_RNA_min_mouse) & (nFeature_RNA >= nFeature_RNA_min_mouse)) %>%
    filter((mitoRatio < mitoRatio_max) & (nCount_RNA <= nCount_RNA_max) & (nFeature_RNA <= nFeature_RNA_max))
  print(nrow(metadata_filtered_df1))
  
  metadata_filtered_df2 <- metadata_filtered_df1 %>%
    filter(!predicted_doublet)
  print(nrow(metadata_filtered_df2))

  ## merge intp the super table
  metadata_seuratfiltered_df <- rbind(metadata_seuratfiltered_df, metadata_filtered_df2)
}

# write unfiltered meta data ----------------------------------------------
file2write <- paste0(dir_out, "CellRanger_Filtered.Barcode_Metrics.tsv")
write.table(x = metadata_cellrangerfiltered_df, file = file2write, sep = "\t", quote = F, row.names = F)

# write seurat filtered meta data ----------------------------------------------
file2write <- paste0(dir_out, "Seurat_Filtered.Barcode_Metrics.tsv")
write.table(x = metadata_seuratfiltered_df, file = file2write, sep = "\t", quote = F, row.names = F)

sink()

