# Yige Wu @WashU Apr 2020
## for make a meta data table for each barcode per sample

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/load_pkgs.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies --------------------------------------------
## set directory to the cell ranger outputs
dir_cellranger_outputs <- "./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/Resources/snRNA_Processed_Data/Cell_Ranger/outputs/GRCh38-3.0.0.premrna_and_mm10_premrna/"
## set sample ids
sample_ids <- list.files(path = dir_cellranger_outputs)
sample_ids

# process by sample -------------------------------------------------------
path_outputs <- NULL
for (sample_id in sample_ids) {
  ## Make matrix directory
  dir_matrix <- paste0(dir_cellranger_outputs, sample_id, "/outs/raw_feature_bc_matrix/")
  
  ## Read in matrix
  input <- Seurat::Read10X(data.dir = dir_matrix)
  srat <- Seurat::CreateSeuratObject(counts = input,
                             min.features = 100, project = sample_id)
  
  ## Compute number of genes detected per UMI
  ### this metric with give us an idea of the complexity of our dataset (more genes detected per UMI, more complex our data)
  srat$log10GenesPerUMI <- log10(srat$nFeature_RNA) / log10(srat$nCount_RNA)
  
  ## Compute percent mito ratio
  srat$mitoRatio <- Seurat::PercentageFeatureSet(object = srat, pattern = "MT-")
  srat$mitoRatio <- srat@meta.data$mitoRatio / 100
  
  # Create metadata dataframe
  metadata_df <- srat@meta.data
  metadata_df$barcode <- rownames(metadata_df)
  
  ## input gem classification file
  path_gem_class <- paste0(dir_cellranger_outputs, sample_id, "/outs/analysis/gem_classification.csv")
  gem_class_df <- fread(input = path_gem_class, data.table = F)
  
  ## change format
  gem_class_df <- gem_class_df %>%
    rename(barcode_raw = barcode) %>%
    mutate(barcode = str_split_fixed(string = barcode_raw, pattern = "-", n = 2)[,1])
  
  ## add gem class info
  metadata_df <- merge(metadata_df, gem_class_df, by = c("barcode"), all.x = T)
  
  ## write output
  file2write <- paste0(dir_out, sample_id, "_quality_metrics_per_barcode.", run_id, ".tsv")
  write.table(x = metadata_df, file = file2write, quote = F, sep = "\t", row.names = F)
  
  ## store output paths
  path_outputs <- c(path_outputs, file2write)
}

# make and write path to the quality metrics ------------------------------
path_outputs_df <- data.frame(sample_id = sample_ids, path_output = path_outputs)
file2write <- paste0(dir_out, "path_quality_metrics_per_barcode.", run_id, ".tsv")
write.table(x = path_outputs_df, file = file2write, quote = F, sep = "\t", row.names = F)
