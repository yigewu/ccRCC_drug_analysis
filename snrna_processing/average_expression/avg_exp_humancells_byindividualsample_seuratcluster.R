# Yige Wu @WashU Apr 2020

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
## input the integrated data
path_rds <- "./Resources/Analysis_Results/snrna_processing/integration/run_humancells_8sample_integration_withanchor_on_katmai/20210208.v1/Humancells_8sample_integration.withanchor.20210208.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-to-tumorsubcluster table
barcode2tumorsubcluster_df <- fread(input = "./Resources/Analysis_Results/snrna_processing/fetch_data/fetch_data_individual_sample_humancells_on_katmai/20210224.v1/HumanCells.Individual_sample.umap_data.20210224.v1.tsv", data.table = F)
# nrow(barcode2tumorsubcluster_df)
cat("finish reading the barcode-to-tumorsubcluster table!\n")
# spcify assay
assay_process <- "SCT"
slot_process <- "scale.data"
cat(paste0("Assay: ", assay_process, "\n"))

# modify srat object ------------------------------------------------------
## add unique barcode id to the seurat meta data
### extract the meta data 
metadata_df <- srat@meta.data
metadata_df$barcode_31integration <- rownames(metadata_df)
metadata_df <- metadata_df %>%
  mutate(barcode_individual = str_split_fixed(string = barcode_31integration, pattern = "_", n = 2)[,1]) %>%
  mutate(aliquot_barcode = paste0(orig.ident, "_", barcode_individual))
### modify the tumor barcode table
barcode2tumorsubcluster_df <- barcode2tumorsubcluster_df %>%
  mutate(aliquot_barcode = paste0(orig.ident, "_", barcode)) %>%
  mutate(Name_TumorSubcluster = paste0(orig.ident, "_C", seurat_clusters))
### add new column for tumor cluster name
tmp <- mapvalues(x = metadata_df$aliquot_barcode, from = barcode2tumorsubcluster_df$aliquot_barcode, to = as.vector(barcode2tumorsubcluster_df$Name_TumorSubcluster))
table(tmp)
## change ident
srat@meta.data$Name_TumorSubcluster <- tmp
Idents(srat) <- "Name_TumorSubcluster"

## run average expression
aliquot.averages <- AverageExpression(srat, assays = assay_process, slot = slot_process)
print("Finish running AverageExpression!\n")
cat("###########################################\n")

## write output
file2write <- paste0(dir_out, "AverageExpression.", "ByIndividualSampleSeuratCluster.", run_id, ".tsv")
write.table(aliquot.averages, file = file2write, quote = F, sep = "\t", row.names = T)
cat("Finished saving the output\n")
cat("###########################################\n")
sink()
