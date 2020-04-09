# Yige Wu @WashU Apr 2020
## run SCTransform by for loop

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/load_pkgs.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/functions.R")
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
path_srat_obj_df <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/Resources/Analysis_Results/snrna_processing/filtering/make_seuratfiltered_srat_obj/20200402.v1/path_seurat_filtered_RDS.20200402.v1.tsv", data.table = F)
# cell cycle genes --------------------------------------------------------
cell_cycle_genes <- fread(input = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/Resources/Analysis_Results/dependencies/make_cellcycle_human_mouse_genes/20200408.v1/cell_cycle_human_mouse_genes.20200408.v1.tsv", data.table = F)
cell_cycle_genes <- cell_cycle_genes %>%
  mutate(feature_name = ifelse(species == "human", paste0("GRCh38-3.0.0.premrna-", gene_name), paste0("mm10-premrna---------", gene_name)))
g2m_feature_names <- cell_cycle_genes$feature_name[cell_cycle_genes$phase == "G2/M"]
s_feature_names <- cell_cycle_genes$feature_name[cell_cycle_genes$phase == "S"]

# process each sample -----------------------------------------------------
path_outputs <- NULL
for (sample_id in path_srat_obj_df$sample_id) {
  ## input filtered srat object
  path_srat <- path_srat_obj_df$path_output[path_srat_obj_df$sample_id == sample_id]
  srat <- readRDS(file = path_srat)
  srat <- NormalizeData(srat, verbose = TRUE)
  srat <- CellCycleScoring(srat, g2m.features=g2m_feature_names, s.features=s_feature_names)
  srat <- SCTransform(srat, vars.to.regress = c("mitoRatio", 'nFeature_RNA', "nCount_RNA", 'S.Score', 'G2M.Score'))
  
  ## save outputs
  file2write <- paste0(dir_out, sample_id, ".seruat_sctransform.", run_id, ".RDS")
  saveRDS(object = srat, file = file2write, compress = T)
  
  ## store path to the outputs
  path_outputs <- c(path_outputs, file2write)
}
# make and write path to the quality metrics ------------------------------
path_outputs_df <- data.frame(sample_id = path_srat_obj_df$sample_id, path_output = path_outputs)
file2write <- paste0(dir_out, "path_seurat_sctransform_RDS.", run_id, ".tsv")
write.table(x = path_outputs_df, file = file2write, quote = F, sep = "\t", row.names = F)

