# Yige Wu @WashU Jan 2021
## filter the Cell Ranger filtered barcode by the same threshold
## create the seurat object

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
metrics_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/filtering/make_seuratfiltered_srat_obj/20210125.v1/worklog.Seurat_Filtered.Barcode_Metrics.tsv")

# summarize ---------------------------------------------------------------
summary_df <- metrics_df %>%
  group_by(orig.ident, call) %>%
  summarise(count_cells = n(), median_nCount_RNA = median(nCount_RNA), median_nFeature_RNA = median(nFeature_RNA))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Seurat_Filtered.", "Quality_Metric_Summary.", run_id, ".tsv")
write.table(file = file2write, x = summary_df, quote = F, sep = "\t", row.names = T)