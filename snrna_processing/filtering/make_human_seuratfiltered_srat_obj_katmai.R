# Yige Wu @WashU Jan 2021
## filter the Cell Ranger filtered barcode by the same threshold
## create the seurat object

# set up libraries and output directory -----------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set time stamp for log file
timestamp <- paste0(run_id, ".", format(Sys.time(), "%H%M%S"))
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
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
library(Seurat)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)
dir_out_elbowplot <- paste0(dir_out, "elbowplot", "/")
dir.create(dir_out_elbowplot)
dir_out_umap <- paste0(dir_out, "umap", "/")
dir.create(dir_out_umap)
## set log file
sink(file = paste0(dir_out, "Log.", timestamp, ".txt"))

# input dependencies ------------------------------------------------------
## set directory to the cell ranger outputs
# dir_cellranger_outputs <- "./Resources/snRNA_Processed_Data/Cell_Ranger/cellranger-5.0.1_Ref-2020-A_GRCh38/"
dir_cellranger_outputs <- "./Resources/snRNA_Processed_Data/Cell_Ranger/outputs/cellranger-5.0.1_Ref-2020-A_GRCh38/"
## set sample ids
ids_sample <- list.files(path = dir_cellranger_outputs)
ids_sample <- ids_sample[grepl(pattern = 'RESL', x = ids_sample)]
ids_sample
## input filtered barcode table
barcode2species_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/map_barcode/map_barcode2celltype_with_8sample_integration/20210129.v1/Barcode2CellType_UMAP.20210129.v1.tsv")
## input cell cycle genes
cell_cycle_genes <- fread(input = "./Resources/Analysis_Results/dependencies/make_cellcycle_human_mouse_genes/20200408.v1/cell_cycle_human_mouse_genes.20200408.v1.tsv", data.table = F)
cell_cycle_genes <- cell_cycle_genes %>%
  mutate(feature_name = ifelse(species == "human", paste0("GRCh38-", gene_name), paste0("mm10---", gene_name)))
g2m_feature_names <- cell_cycle_genes$feature_name[cell_cycle_genes$phase == "G2/M"]; print(g2m_feature_names)
s_feature_names <- cell_cycle_genes$feature_name[cell_cycle_genes$phase == "S"]; print(s_feature_names)

# process each sample -----------------------------------------------------
path_outputs <- NULL
for (id_sample in ids_sample) {
  print(id_sample)
  ## get the barcodes to keep
  barcodes_keep_df <- barcode2species_df %>%
    filter(Id_Sample == id_sample) %>%
    filter(Species_Cell == "Human") %>%
    mutate(Barcode_Raw = str_split_fixed(string = Barcode_Integrated, pattern = "_", n = 2)[,1])
  barcodes_keep <- barcodes_keep_df$Barcode_Raw
  ## Make matrix directory
  dir_matrix <- paste0(dir_cellranger_outputs, id_sample, "/", id_sample, "/outs/filtered_feature_bc_matrix/")
  
  ## Read in matrix
  input <- Seurat::Read10X(data.dir = dir_matrix); print(dim(input))
  
  ## subset the matrix
  input_subset <- input[, barcodes_keep]; print(dim(input_subset))
  rm(input)
  
  ## make seurat object
  srat <- Seurat::CreateSeuratObject(counts = input_subset, project = id_sample)
  rm(input_subset)
  
  ## normalize
  srat <- NormalizeData(srat, verbose = TRUE)
  
  ## add cell cycle score
  print(tail(srat@assays$RNA@data@Dimnames[[1]]))
  print(head(srat@assays$RNA@data@Dimnames[[1]]))
  srat <- CellCycleScoring(object = srat, g2m.features=g2m_feature_names, s.features=s_feature_names)
  
  ## SCTransform
  srat <- SCTransform(srat, vars.to.regress = c("mitoRatio", 'nFeature_RNA', "nCount_RNA", 'S.Score', 'G2M.Score'))
  
  # Run PCA
  srat <- RunPCA(object = srat, npcs = num_pcs)
  
  # Plot ans save the elbow plot, see if the number of PCs are reasonable
  p <- ElbowPlot(object = srat, ndims = num_pcs)
  file2write <- paste0(dir_out_elbowplot, id_sample, "_pc_elbowplot", ".png")
  png(filename = file2write, width = 800, height = 800, res = 150)
  print(p)
  dev.off()
  
  # Run UMAP
  srat <- RunUMAP(srat, dims = 1:num_pcs, reduction = "pca")
  
  # Determine the K-nearest neighbor graph
  srat <- FindNeighbors(object = srat, dims = 1:num_pcs)
  
  # Determine the clusters for various resolutions                                
  srat <- FindClusters(object = srat, resolution = findclusters_res)
  
  # Plot ans save the UMAP plot
  p <- DimPlot(srat,
               reduction = "umap",
               label = TRUE,
               label.size = 6)
  file2write <- paste0(dir_out_umap, id_sample, "_", paste0("SCT_snn_res.", findclusters_res), ".umap.", ".png")
  png(filename = file2write, width = 1000, height = 800, res = 150)
  print(p)
  dev.off()
  
  ## save outputs
  file2write <- paste0(dir_out, id_sample, ".human", ".seruat_filtered.", run_id, ".RDS")
  saveRDS(object = srat, file = file2write, compress = T)
  
  ## store path to the outputs
  path_outputs <- c(path_outputs, file2write)
}
# make and write path to the quality metrics ------------------------------
path_outputs_df <- data.frame(id_sample = ids_sample, 
                              path_output_katmai = path_outputs, 
                              path_output_box = gsub(x = path_outputs, pattern = dir_base, replacement = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"), 
                              path_output_relative = gsub(x = path_outputs, pattern = dir_base, replacement = "./"))
file2write <- paste0(dir_out, "Path_to_Seurat_Objects.HumanCellsss.Filtered.", run_id, ".tsv")
write.table(x = path_outputs_df, file = file2write, quote = F, sep = "\t", row.names = F)
sink()

