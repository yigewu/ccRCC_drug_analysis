# Yige Wu @WashU Feb 2021

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
## input mouse cell type
mouse_bc2ct_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/map_barcode/map_barcode2celltype_mousecells_integration_withanchor/20210225.v1/MouseCells.Barcode2CellType.20210225.v1.tsv")
## input human cell type
human_bc2ct_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/map_barcode/map_barcode2celltype_humancells_integration_withanchor/20210225.v1/HumanCells.Barcode2CellType.20210225.v1.tsv")

# combine -----------------------------------------------------------------
mouse_bc2ct_df2merge <- mouse_bc2ct_df %>%
  rename(Cell_Group4 = Cell_Group3) %>%
  rename(Cell_Group3 = Cell_Group2) %>%
  mutate(Cell_Species = "Mouse")
human_bc2ct_df2merge <- human_bc2ct_df %>%
  rename(Cell_Group3 = Cell_Group) %>%
  mutate(Cell_Group4 = Cell_Group3) %>%
  mutate(Cell_Species = "Human")
merged_bc2ct_df <- rbind(mouse_bc2ct_df2merge, human_bc2ct_df2merge)
merged_bc2ct_df <- merged_bc2ct_df %>%
  select(Cell_Species, Id_Sample, 
         Cell_Type.Short, Cell_Type.Detailed, Cell_Group3, Cell_Group4,
         Barcode_Integrated, UMAP_1, UMAP_2, Id_Seurat_Cluster, Treatment_Group)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Human_and_Mouse.", "Barcode2CellType.", run_id, ".", "tsv")
write.table(x = merged_bc2ct_df, file = file2write, sep = "\t", row.names = F, quote = F)
