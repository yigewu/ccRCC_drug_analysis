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
## input gene order file
gene_order_df = read.delim(file = "./Resources/snRNA_Processed_Data/InferCNV/inputs/gene_order_file/gencode_v21_gen_pos.MatchCellRangerFeatures.NoDuplicates.20191005.v1.txt", header = FALSE,stringsAsFactors = FALSE)

# process by sample -------------------------------------------------------
for (id_sample_tmp in path_rds_df$id_sample) {
  ## get the path for RDS file
  path_rds <-path_rds_df$path_output_relative[path_rds_df$id_sample == id_sample_tmp]
  ## read RDS file and store in the list
  srat <- readRDS(file = path_rds)
  
  ## get raw read count matrix
  raw_exp_mat <- srat@assays$RNA@counts
  dim(raw_exp_mat)
  
  ## remove gene not mapped in the gene position file
  missing_sym  <- rownames(raw_exp_mat)[!(rownames(raw_exp_mat) %in% gene_order_df$V1)]
  missing_sym
  
  ## only missing ~750 genes, romove them
  ## only keep the barcodes in the annotation files
  clean_exp_mat <- raw_exp_mat[!(rownames(raw_exp_mat) %in% missing_sym), ]
  rm(raw_exp_mat)
  clean_exp_mat <- as.matrix(clean_exp_mat)
  
  path_out <- paste0(dir_out, id_sample_tmp, ".Filtered.Count_Matrix.tsv")
  write.table(x = clean_exp_mat, file = path_out, quote = F, row.names = T, col.names = T, sep = "\t")
}

