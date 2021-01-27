# Yige Wu @WashU Jan 2021
## filter the Cell Ranger filtered barcode by the same threshold
## create the seurat object

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
library(Seurat)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
## set time stamp for log file
timestamp <- paste0(run_id, ".", format(Sys.time(), "%H%M%S"))
## set log file
sink(file = paste0(dir_out, "Log.", timestamp, ".txt"))
## create sub-output directory
dir_out_elbowplot <- paste0(dir_out, "elbowplot", "/")
dir.create(dir_out_elbowplot)
dir_out_umap <- paste0(dir_out, "umap", "/")
dir.create(dir_out_umap)

# input dependencies ------------------------------------------------------
## set directory to the cell ranger outputs
dir_cellranger_outputs <- "./Resources/snRNA_Processed_Data/Cell_Ranger/outputs/cellranger-5.0.1_Ref-2020-A/"
## set sample ids
ids_sample <- list.files(path = dir_cellranger_outputs)
ids_sample <- ids_sample[grepl(pattern = 'RESL', x = ids_sample)]
# ids_sample
## input paths to the scrublet outputs
scrublet_out_paths_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/other/write_paths_to_scrublet_output_tables/20210127.v1/Paths_to_Scrulbet_Output_Tables.20210127.v1.tsv")
## input scrublet 
## set filtering threshold
mitoRatio_max <- 0.1
nCount_RNA_min_human <- 1000
nCount_RNA_min_mouse <- 500
nCount_RNA_max <- 80000
nFeature_RNA_min_human <- 200
nFeature_RNA_min_mouse <- 100
nFeature_RNA_max <- 10000
## input cell cycle genes
cell_cycle_genes <- fread(input = "./Resources/Analysis_Results/dependencies/make_cellcycle_human_mouse_genes/20200408.v1/cell_cycle_human_mouse_genes.20200408.v1.tsv", data.table = F)
cell_cycle_genes <- cell_cycle_genes %>%
  mutate(feature_name = ifelse(species == "human", paste0("GRCh38-", gene_name), paste0("mm10---", gene_name)))
g2m_feature_names <- cell_cycle_genes$feature_name[cell_cycle_genes$phase == "G2/M"]; print(g2m_feature_names)
s_feature_names <- cell_cycle_genes$feature_name[cell_cycle_genes$phase == "S"]; print(s_feature_names)

# process each sample -----------------------------------------------------
metadata_cellrangerfiltered_df <- NULL
path_outputs <- NULL
metadata_seuratfiltered_df <- NULL
for (id_sample in ids_sample) {
  print(id_sample)
  
  ## Make matrix directory
  dir_matrix <- paste0(dir_cellranger_outputs, id_sample, "/", id_sample, "/outs/filtered_feature_bc_matrix/")
  
  ## Read in matrix
  input <- Seurat::Read10X(data.dir = dir_matrix)
  srat <- Seurat::CreateSeuratObject(counts = input,
                                     project = id_sample)
  rm(input)
  
  ## Compute percent mito ratio
  srat$mitoRatio <- Seurat::PercentageFeatureSet(object = srat, pattern = "MT-")
  srat$mitoRatio <- srat@meta.data$mitoRatio / 100
  
  ## extract meta data
  metadata_df <- srat@meta.data
  metadata_df$barcode <- rownames(metadata_df)
  print(head(metadata_df))
  
  ## add barcode classification
  class_df <- fread(data.table = F, input = paste0(dir_cellranger_outputs, id_sample, "/", id_sample, "/outs/analysis/gem_classification.csv"))
  print(head(class_df))
  
  metadata_df <- merge(x = metadata_df, y = class_df, by.x = c("barcode"), by.y = c("barcode"), all.x = T)
  print(head(metadata_df))
  
  ## add scrublet information if any
  ### input scrublet output
  if (id_sample %in% scrublet_out_paths_df$Sample_id) {
    path_scrublet <- scrublet_out_paths_df$Path_relative[scrublet_out_paths_df$Sample_id == id_sample]
    scrublet_out_df <- fread(data.table = F, input = path_scrublet)
    metadata_df <- merge(x = metadata_df, y = scrublet_out_df, by.x = c("barcode"), by.y = c("Barcode"), all.x = T)
  } else {
    metadata_df <- metadata_df %>%
      mutate(doublet_score = NA) %>%
      mutate(predicted_doublet = FALSE)
  }
  print(head(metadata_df))
  ## merge into the super table
  metadata_cellrangerfiltered_df <- rbind(metadata_cellrangerfiltered_df, metadata_df)
  
  ## filter by quality metrics
  metadata_filtered_df1 <- metadata_df %>%
    filter(((call == "mm10") & (nCount_RNA >= nCount_RNA_min_mouse) & (nFeature_RNA >= nFeature_RNA_min_mouse)) | ((call == "GRCh38") & (nCount_RNA >= nCount_RNA_min_human) & (nFeature_RNA >= nFeature_RNA_min_human))) %>%
    filter((mitoRatio < mitoRatio_max) & (nCount_RNA <= nCount_RNA_max) & (nFeature_RNA <= nFeature_RNA_max))
  print(nrow(metadata_filtered_df1))
  
  metadata_filtered_df2 <- metadata_filtered_df1 %>%
    filter(!predicted_doublet)
  print(nrow(metadata_filtered_df2))
  ## subset seurat object
  filtered_srat <- subset(x = srat, cells = metadata_filtered_df2$barcode)
  rm(srat)
  
  ## add meta data
  rownames(metadata_filtered_df2) <- metadata_filtered_df2$barcode
  filtered_srat@meta.data <- metadata_filtered_df2
  
  ## normalize
  filtered_srat <- NormalizeData(filtered_srat, verbose = TRUE)
  
  ## add cell cycle score
  print(tail(filtered_srat@assays$RNA@data@Dimnames[[1]]))
  print(head(filtered_srat@assays$RNA@data@Dimnames[[1]]))
  filtered_srat <- CellCycleScoring(object = filtered_srat, g2m.features=g2m_feature_names, s.features=s_feature_names)
  
  ## SCTransform
  filtered_srat <- SCTransform(filtered_srat, vars.to.regress = c("mitoRatio", 'nFeature_RNA', "nCount_RNA", 'S.Score', 'G2M.Score'))
  
  # Run PCA
  filtered_srat <- RunPCA(object = filtered_srat, npcs = num_pcs)
  
  # Plot ans save the elbow plot, see if the number of PCs are reasonable
  p <- ElbowPlot(object = filtered_srat, ndims = num_pcs)
  file2write <- paste0(dir_out_elbowplot, id_sample, "_pc_elbowplot", ".png")
  png(filename = file2write, width = 800, height = 800, res = 150)
  print(p)
  dev.off()
  
  # Run UMAP
  filtered_srat <- RunUMAP(filtered_srat, dims = 1:num_pcs, reduction = "pca")
  
  # Determine the K-nearest neighbor graph
  filtered_srat <- FindNeighbors(object = filtered_srat, dims = 1:num_pcs)
  
  # Determine the clusters for various resolutions                                
  filtered_srat <- FindClusters(object = filtered_srat, resolution = findclusters_res)
  
  # Plot ans save the UMAP plot
  p <- DimPlot(filtered_srat,
               reduction = "umap",
               label = TRUE,
               label.size = 6)
  file2write <- paste0(dir_out_umap, id_sample, "_", paste0("SCT_snn_res.", findclusters_res), ".umap.", ".png")
  png(filename = file2write, width = 1000, height = 800, res = 150)
  print(p)
  dev.off()
  
  ## save outputs
  file2write <- paste0(dir_out, id_sample, ".seruat_filtered.", run_id, ".RDS")
  saveRDS(object = filtered_srat, file = file2write, compress = T)
  
  ## store path to the outputs
  path_outputs <- c(path_outputs, file2write)
  
  ## merge intp the super table
  metadata_seuratfiltered_df <- rbind(metadata_seuratfiltered_df, filtered_srat@meta.data)
}
# make and write path to the quality metrics ------------------------------
path_outputs_df <- data.frame(id_sample = ids_sample, 
                              path_output_katmai = path_outputs, 
                              path_output_box = gsub(x = path_outputs, pattern = dir_base, replacement = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"), 
                              path_output_relative = gsub(x = path_outputs, pattern = dir_base, replacement = "./"))
file2write <- paste0(dir_out, "path_seurat_filtered_RDS.", run_id, ".tsv")
write.table(x = path_outputs_df, file = file2write, quote = F, sep = "\t", row.names = F)

# write unfiltered meta data ----------------------------------------------
file2write <- paste0(dir_out, "CellRanger_Filtered.Barcode_Metrics.tsv")
write.table(x = metadata_cellrangerfiltered_df, file = file2write, sep = "\t", quote = F, row.names = F)

# write seurat filtered meta data ----------------------------------------------
file2write <- paste0(dir_out, "Seurat_Filtered.Barcode_Metrics.tsv")
write.table(x = metadata_seuratfiltered_df, file = file2write, sep = "\t", quote = F, row.names = F)

sink()

