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
sigScores <- readRDS("./Resources/Analysis_Results/snrna_processing/signature_scores/run_vision/run_vision_on_humancells_8sample_SCT_counts/20221216.v1/humancells.8sampleintegrated.Vision.20221216.v1.RDS")
print("Finish reading the sigScores matrix!\n")
## input the barcode to new cluster
barcode2cluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_byres_humancells_8sample_integration_on_katmai/20220905.v1/HumanCells.8sample.Metadata.ByResolution.20220905.v1.tsv")

# prepare data to plot -----------------------------------------------------------------
barcode2cluster_df <- barcode2cluster_df %>%
  mutate(clusterid_test = integrated_snn_res.0.5) %>%
  mutate(model_id = str_split_fixed(string = orig.ident, pattern = "\\-", n = 3)[,1]) %>%
  mutate(group = paste0(orig.ident, "-", clusterid_test))
clusters <- unique(barcode2cluster_df$clusterid_test)
genesets_test <- colnames(sigScores)
barcodes_all <- barcode2cluster_df$barcode
groups_all <- unique(barcode2cluster_df$group)
registerDoParallel(cores = 20)

print("Start foreach!\n")
start_time <- Sys.time()
result_list<-foreach(g=genesets_test) %dopar% {
  sigscores_df <- barcode2cluster_df
  sigscores_df$value <- sigScores[barcodes_all, g]
  sigscores_df2 <- sigscores_df %>%
    group_by(group) %>%
    summarise(value = median(value))
  sigscores_df2$gene_set <- g
  rownames(sigscores_df2) <- sigscores_df2$group
  return(sigscores_df2)
}
print("Finish foreach!\n")
results_df <- do.call(rbind.data.frame, result_list)

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "median_scores.res05.bycluster.bysample.humancells.8sampleintegrated.", run_id,".tsv")
write.table(x = results_df, file = file2write, quote = F, sep = "\t", row.names = F)



