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

# input dependencies ------------------------------------------------------
## input median(?) signature scores per cluster
results_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/signature_scores/compare_scores_across_clusters/wilcox_eachcluster_vs_rest_bygeneset_byres05/20220908.v1/wilcox.eachcluster_vs_rest.res1.humancells.8sampleintegrated.20220908.v1.tsv")

# make matrix data for heatmap body color-----------------------------------------------------------------
## extract the data for the matrix
plotdata_df <- results_df %>%
  mutate(x_plot = gene_set) %>%
  mutate(y_plot = cluster) %>%
  mutate(value = median_sigScore)
plotdata_wide_df <- dcast(data = plotdata_df, formula = x_plot ~ y_plot, value.var = "value")
colnames(plotdata_wide_df)[1] <- "gene_set"
plotdata_mat <- as.matrix(plotdata_wide_df[,-1])
rownames(plotdata_mat) <- plotdata_wide_df$gene_set

# write output ------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
## write output
file2write <- paste0(dir_out, "Median_signature_scores_per_res0.5_cluster.", run_id, ".RDS")
saveRDS(object = plotdata_mat, file = file2write, compress = T)


