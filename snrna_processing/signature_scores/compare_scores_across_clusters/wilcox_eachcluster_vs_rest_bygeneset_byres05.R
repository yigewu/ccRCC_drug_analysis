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
if (!("VISION" %in% installed.packages()[,1])) {
  print(paste0(pkg_name_tmp, "is being installed!"))
  devtools::install_github("YosefLab/VISION", dependencies = T)
}
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
barcode2cluster_df$clusterid_test <- barcode2cluster_df$integrated_snn_res.0.5
clusters <- unique(barcode2cluster_df$clusterid_test)
genesets_test <- colnames(sigScores)
registerDoParallel(cores = 20)

results_df <- NULL
for (cluster1 in clusters) {
  cat(paste0("Start compare, ", cluster1, " vs. ", "the rest", "!\n"))
  
  barcodes_cluster1 <- barcode2cluster_df$barcode[barcode2cluster_df$clusterid_test == cluster1]
  barcodes_cluster2 <- barcode2cluster_df$barcode[barcode2cluster_df$clusterid_test != cluster1]
  
  print("Start foreach!\n")
  start_time <- Sys.time()
  result_list<-foreach(g=genesets_test) %dopar% {
    sigscores_cluster1 <- sigScores[barcodes_cluster1, g]
    sigscores_cluster2 <- sigScores[barcodes_cluster2, g]
    median_cluster1 <- median(sigscores_cluster1)
    median_diff <- median(sigscores_cluster1) - median(sigscores_cluster2)
    log2FC <- log2(median(sigscores_cluster1)/median(sigscores_cluster2))
    stat <- wilcox.test(x = sigscores_cluster1, y = sigscores_cluster2)
    p_val <- stat$p.value
    result_tmp <- list(c(p_val, median_diff, log2FC, median_cluster1))
    return(result_tmp)
  }
  print("Finish foreach!\n")
  
  end_time <- Sys.time()
  end_time - start_time 
  result_tmp_df <- data.frame(matrix(data = unlist(result_list), ncol = 4, byrow = T))
  print("Finish rbind.data.frame result_list!\n")
  print(head(result_tmp_df))
  
  colnames(result_tmp_df) <- c("p_val", "median_diff", "log2FC", "median_sigScore")
  result_tmp_df$fdr <- p.adjust(p = result_tmp_df$p_val, method = "fdr")
  result_tmp_df$gene_set <- genesets_test
  result_tmp_df$cluster <- cluster1
  print("Finish adding id columns!\n")
  
  results_df <- rbind(results_df, result_tmp_df)
}

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "wilcox.eachcluster_vs_rest.res1.humancells.8sampleintegrated.", run_id,".tsv")
write.table(x = results_df, file = file2write, quote = F, sep = "\t", row.names = F)



