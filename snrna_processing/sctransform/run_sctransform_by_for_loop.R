# Yige Wu @WashU Apr 2020
## run SCTransform by for loop

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
library(Seurat)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the path to the seurat-filtered srat objects
path_srat_obj_df <- fread("./Resources/Analysis_Results/snrna_processing/filtering/make_seuratfiltered_srat_obj/20200416.v1/path_seurat_filtered_RDS.20200416.v1.tsv", data.table = F)
# cell cycle genes --------------------------------------------------------
cell_cycle_genes <- fread(input = "./Resources/Analysis_Results/dependencies/make_cellcycle_human_mouse_genes/20200408.v1/cell_cycle_human_mouse_genes.20200408.v1.tsv", data.table = F)
cell_cycle_genes <- cell_cycle_genes %>%
  mutate(feature_name = ifelse(species == "human", paste0("GRCh38-3.0.0.premrna-", gene_name), paste0("mm10-premrna---------", gene_name)))
g2m_feature_names <- cell_cycle_genes$feature_name[cell_cycle_genes$phase == "G2/M"]
s_feature_names <- cell_cycle_genes$feature_name[cell_cycle_genes$phase == "S"]

# process each sample -----------------------------------------------------
path_outputs <- NULL
for (id_sample in path_srat_obj_df$id_sample) {
  ## input filtered srat object
  path_srat <- path_srat_obj_df$path_output_relative[path_srat_obj_df$id_sample == id_sample]
  srat <- readRDS(file = path_srat)
  srat <- NormalizeData(srat, verbose = TRUE)
  srat <- CellCycleScoring(srat, g2m.features=g2m_feature_names, s.features=s_feature_names)
  srat <- SCTransform(srat, vars.to.regress = c("mitoRatio", 'nFeature_RNA', "nCount_RNA", 'S.Score', 'G2M.Score'))
  
  ## save outputs
  file2write <- paste0(dir_out, id_sample, ".seruat_sctransform.", run_id, ".RDS")
  saveRDS(object = srat, file = file2write, compress = T)
  
  ## store path to the outputs
  path_outputs <- c(path_outputs, file2write)
}
# make and write path to the quality metrics ------------------------------
path_outputs_df <- data.frame(id_sample = path_srat_obj_df$id_sample, 
                              path_output_local = path_outputs,
                              path_output_relative = gsub(x = path_outputs, pattern = dir_base, replacement = "./"))
file2write <- paste0(dir_out, "path_seurat_sctransform_RDS.", run_id, ".tsv")
write.table(x = path_outputs_df, file = file2write, quote = F, sep = "\t", row.names = F)

