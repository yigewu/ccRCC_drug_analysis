# Yige Wu @WashU Map 2020

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
## input the barcode2celltype
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/snrna_processing/fetch_data/fetch_data_RESL_8sample_integration_withanchor_on_katmai/20210128.v1/RESL_8sample.umap_data.20210128.v1.tsv", data.table = F)

# count by cluster by call ------------------------------------------------
count_bycluster_bycall_df <- barcode2celltype_df %>%
  select(seurat_clusters, call) %>%
  table() %>%
  data.frame() %>%
  arrange(seurat_clusters)

