# Yige Wu @WashU Apr 2020
## check if cell cycle phase is a significant variation that needs to be regressed out before sctransform

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/load_pkgs.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/functions.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/plotting.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/variables.R")

library(Seurat)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the path to the seurat-filtered srat objects
path_srat_obj_df <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/Resources/Analysis_Results/quality_control/filtering/make_seuratfiltered_srat_obj/20200402.v1/path_seurat_filtered_RDS.20200402.v1.tsv", data.table = F)
# Load cell cycle markers
cell_cycle_genes <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/Resources/Analysis_Results/dependencies/make_cellcycle_human_mouse_genes/20200408.v1/cell_cycle_human_mouse_genes.20200408.v1.tsv", data.table = F)


# add reference name to human and mouse gene names -----------------------------
## seurat_phase@assays$RNA@counts %>% rownames() %>% tail(n = 1000)
## seurat_phase@assays$RNA@counts %>% rownames() %>% head(n = 1000)
cell_cycle_genes <- cell_cycle_genes %>%
  mutate(feature_name = ifelse(species == "human", paste0("GRCh38-3.0.0.premrna-", gene_name), paste0("mm10-premrna---------", gene_name)))
## check if the cell cycle genes are in the feature names
cell_cycle_genes$feature_name[cell_cycle_genes$phase == "S"] %in% rownames(seurat_phase@assays$RNA@counts)
cell_cycle_genes$feature_name[cell_cycle_genes$phase == "G2/M"] %in% rownames(seurat_phase@assays$RNA@counts)

# process each sample -----------------------------------------------------
for (sample_id in path_srat_obj_df$sample_id) {
  ## input filtered srat object
  path_srat <- path_srat_obj_df$path_output[path_srat_obj_df$sample_id == sample_id]
  srat <- readRDS(file = path_srat)
  
  # Normalize the counts
  seurat_phase <- NormalizeData(object = srat)
  # Score cells for cell cycle
  seurat_phase <- CellCycleScoring(seurat_phase,
                                   g2m.features = cell_cycle_genes$feature_name[cell_cycle_genes$phase == "G2/M"], 
                                   s.features = cell_cycle_genes$feature_name[cell_cycle_genes$phase == "S"])
  
  # View cell cycle scores and phases assigned to cells                                 
  View(seurat_phase@meta.data)  
  
  # Identify the most variable genes
  seurat_phase <- FindVariableFeatures(seurat_phase, 
                                       selection.method = "vst",
                                       nfeatures = num_var_features, 
                                       verbose = FALSE)
  
  # Scale the counts
  seurat_phase <- ScaleData(seurat_phase)
  
  # Perform PCA
  seurat_phase <- RunPCA(seurat_phase, npcs = num_pcs)
  
  # Plot the PCA colored by cell cycle phase
  p <- DimPlot(seurat_phase,
          reduction = "pca",
          group.by= "call",
          split.by = "Phase")
  p
  ## save output
  file2write <- paste0(dir_out, sample_id, "top_pcs_by_cellcyclephase.", run_id, ".png")
  png(filename = file2write, width = 1600, height = 500)
  print(p)
  dev.off()
}
