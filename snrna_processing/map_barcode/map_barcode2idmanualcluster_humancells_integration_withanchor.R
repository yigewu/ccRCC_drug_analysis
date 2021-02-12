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
fetcheddata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/fetch_data/fetch_data_humancells_8sample_integration_withanchor_on_katmai/20210209.v1/HumanCells_8sample.umap_data.20210209.v1.tsv")
## input cell type assignment cluster 2 cell type
seuratcluster2manualcluster_df <- readxl::read_excel(path = "./Resources/snRNA_Processed_Data/Manual_Grouping/HumanCells_Integration_Manual_Clustering.xlsx")

# map barcode2celltype ----------------------------------------------------
fetcheddata_df <- merge(fetcheddata_df, 
                        seuratcluster2manualcluster_df,
                        by.x = c("seurat_clusters"),
                        by.y = c("Id_Seurat_Cluster"), all.x = T)
fetcheddata_df <- fetcheddata_df %>%
  dplyr::rename(Id_Seurat_Cluster = seurat_clusters) %>%
  dplyr::rename(Barcode_Integrated = V1) %>%
  dplyr::rename(Id_Sample = orig.ident) %>%
  dplyr::select(-barcode) %>%
  dplyr::mutate(Treatment_Group = str_split_fixed(string = Id_Sample, pattern = "-", n = 3)[,3]) %>%
  dplyr::mutate(Treatment_Group = gsub(pattern = "2", replacement = "", x = Treatment_Group))
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "barcode2idmanualcluster_umapdata.", run_id, ".tsv")
write.table(x = fetcheddata_df, file = file2write, quote = F, row.names = F, sep = "\t")

