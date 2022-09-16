# Yige Wu @WashU Sep 2022
## before run script, run conda activate vision

# set up libraries and output directory -----------------------------------
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
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "doParallel"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_drug_analysis/functions.R")
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input  ------------------------------------------------------
## input signature scores by barcode
sigScores <- readRDS("./Resources/Analysis_Results/snrna_processing/signature_scores/run_vision/run_vision_on_humancells_8sample_integrated/20220907.v1/humancells.8sampleintegrated.Vision.sigScores.20220907.v1.RDS")
print("Finish reading the sigScores matrix!\n")
## input the barcode to new cluster
barcode2cluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_byres_humancells_8sample_integration_on_katmai/20220905.v1/HumanCells.8sample.Metadata.ByResolution.20220905.v1.tsv")

# prepare data to plot -----------------------------------------------------------------
barcode2cluster_df <- barcode2cluster_df %>%
  mutate(clusterid_test = integrated_snn_res.0.5) %>%
  mutate(model_id = str_split_fixed(string = orig.ident, pattern = "\\-", n = 3)[,1])
clusters <- unique(barcode2cluster_df$clusterid_test)
genesets_test <- colnames(sigScores)
registerDoParallel(cores = 20)

results_df <- NULL
for (cluster_tmp in clusters) {
  barcodes_group1 <- barcode2cluster_df$barcode[barcode2cluster_df$model_id == "RESL5E" & barcode2cluster_df$clusterid_test == cluster_tmp]
  barcodes_group2 <- barcode2cluster_df$barcode[barcode2cluster_df$model_id == "RESL10F" & barcode2cluster_df$clusterid_test == cluster_tmp]
  
  print("Start foreach!\n")
  start_time <- Sys.time()
  result_list<-foreach(g=genesets_test) %dopar% {
    sigscores_group1 <- sigScores[barcodes_group1, g]
    sigscores_group2 <- sigScores[barcodes_group2, g]
    median_group1 <- median(sigscores_group1)
    median_group2 <- median(sigscores_group2)
    
    result_tmp <- list(c(g, median_group1, median_group2))
    return(result_tmp)
  }
  print("Finish foreach!\n")
  result_tmp_df <- data.frame(matrix(data = unlist(result_list), ncol = 3, byrow = T))
  result_tmp_df$clusterid_test <- cluster_tmp
  results_df <- rbind(results_df, result_tmp_df)
}

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "median_scores.res05.humancells.8sampleintegrated.", run_id,".tsv")
write.table(x = results_df, file = file2write, quote = F, sep = "\t", row.names = F)



