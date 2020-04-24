# Yige Wu @WashU Apr 2020

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
## input fetched data
fetcheddata_df <- fread(input = "./Resources/Analysis_Results/snrna_processing/fetch_data/fetch_data_RESL5_4sample_integration_withanchor_on_katmai/20200423.v1/RESL5_4sample_integration.withanchor.20200417.v1.fetched_data.20200423.v1.tsv", data.table = F)
fetcheddata_df$Id_Model <- "RESL5"
## input cell type assignment cluster 2 cell type
cluster2celltype_df <- readxl::read_excel(path = "./Resources/snRNA_Processed_Data/Cell_Type_Assignment/4sample_Integration/4sample_Integration_Cluster2CellType.20200423.v1.xlsx")

# map barcode2celltype ----------------------------------------------------
fetcheddata_df <- merge(fetcheddata_df, cluster2celltype_df, 
                        by.x = c("Id_Model", "seurat_clusters"),
                        by.y = c("Id_Model", "Cluster"), all.x = T)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "barcode2celltype_umapdata.", run_id, ".tsv")
write.table(x = fetcheddata_df, file = file2write, quote = F, row.names = F, sep = "\t")

