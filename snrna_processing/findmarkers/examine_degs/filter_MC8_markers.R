# Yige Wu @WashU 2023

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_unique_markers_for_selected_clusters_katmai/20230411.v1/MC8.logfcthreshold.0.25.minpct.0.1.mindiffpct.0.tsv")

# count DEG across samples ------------------------------------------------
deg_count_df = deg_df %>%
  filter(p_val_adj < 0.05) %>%
  select(gene_symbol) %>%
  table() %>%
  as.data.frame()

