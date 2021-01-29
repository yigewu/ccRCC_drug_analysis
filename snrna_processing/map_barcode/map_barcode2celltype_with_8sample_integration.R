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
fetcheddata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/fetch_data/fetch_data_RESL_8sample_integration_withanchor_on_katmai/20210128.v1/RESL_8sample.umap_data.20210128.v1.tsv")
## input cell type assignment cluster 2 cell type
cluster2celltype_df <- readxl::read_excel(path = "./Resources/snRNA_Processed_Data/Cell_Type_Assignment/8sample_Integration/8sample_Integration_Cluster2CellType.20210129.v1.xlsx")

# map barcode2celltype ----------------------------------------------------
fetcheddata_df <- merge(fetcheddata_df %>%
                          select(-V1), 
                        cluster2celltype_df %>%
                          select(-Comment), 
                        by.x = c("seurat_clusters", "call"),
                        by.y = c("Cluster", "Call"), all.x = T)
fetcheddata_df <- fetcheddata_df %>%
  rename(Id_Seurat_Cluster = seurat_clusters) %>%
  rename(Barcode_Integrated = barcode) %>%
  rename(Id_Sample = orig.ident)
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Barcode2CellType_UMAP.", run_id, ".tsv")
write.table(x = fetcheddata_df, file = file2write, quote = F, row.names = F, sep = "\t")

