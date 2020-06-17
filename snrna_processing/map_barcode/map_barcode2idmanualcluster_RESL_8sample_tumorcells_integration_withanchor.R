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
fetcheddata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/clustering/cluster_RESL_tumorcells_integration_withanchor_on_katmai/20200507.v1/umap_data.RESL.Tumor_cells.integration.withanchor.20200507.v1.PC40.Res0.5.tsv")
## input cell type assignment cluster 2 cell type
seuratcluster2manualcluster_df <- readxl::read_excel(path = "./Resources/snRNA_Processed_Data/Manual_Grouping/manual_cluster_grouping.RESL_8sample_tumorcells_integration_withanchor.20200508.v1.xlsx")

# map barcode2celltype ----------------------------------------------------
fetcheddata_df <- merge(fetcheddata_df, 
                        seuratcluster2manualcluster_df,
                        by.x = c("seurat_clusters"),
                        by.y = c("Id_Seurat_Cluster"), all.x = T)
fetcheddata_df <- fetcheddata_df %>%
  rename(Id_Seurat_Cluster = seurat_clusters) %>%
  rename(Barcode_Integrated = V1) %>%
  rename(Id_Sample = orig.ident)
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "barcode2idmanualcluster_umapdata.", run_id, ".tsv")
write.table(x = fetcheddata_df, file = file2write, quote = F, row.names = F, sep = "\t")

