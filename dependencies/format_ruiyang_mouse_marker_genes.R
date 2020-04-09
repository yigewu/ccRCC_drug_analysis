# Yige Wu @WashU Apr 2020
## format mouse kidney markers from Ruiyang to Yige's format


# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory

dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input Ruiyang's marker gene lsit
genes_celltypemarker_df <- fread(input = "./Resources/Gene_Lists/Cell_Type_Marker_Genes/KidneyCellMarkers.txt", data.table = F)

# examine gene list -------------------------------------------------------
## take out unwanted markers
genes_celltypemarker_df <- genes_celltypemarker_df %>%
  filter(!(cell_type %in% c("PKD", "QC")))
genes_celltypemarker_df <- genes_celltypemarker_df %>%
  rename(Cell_Type_Short = cell_type) %>%
  rename(Gene_Symbol = gene_symbol)
genes_celltypemarker_df %>%
  select(Cell_Type_Short) %>%
  table()

# add cell type groups -----------------------------------------------------
genes_celltypemarker_df$Cell_Type_Group <- "Other"
genes_celltypemarker_df$Cell_Type_Group[genes_celltypemarker_df$Cell_Type_Short %in% c("Urothelium")] <- "Urothelium"
genes_celltypemarker_df$Cell_Type_Group[genes_celltypemarker_df$Cell_Type_Short %in% c("Distalconvolutedtubules", "Intercalatedcells", "LoopofHenle", "Podocytes", "PrincipleCells", "Proximaltubules")] <- "Neprhon_Epithelium"
genes_celltypemarker_df$Cell_Type_Group[genes_celltypemarker_df$Cell_Type_Short %in% c("Endothelial", "Fibroblasts", "Myofibroblasts")] <- "Stroma"
genes_celltypemarker_df$Cell_Type_Group[genes_celltypemarker_df$Cell_Type_Short %in% c("Macrophages")] <- "Immune"
genes_celltypemarker_df %>%
  select(Cell_Type_Group) %>%
  table()
## add cell type 1
genes_celltypemarker_df$Cell_Type1 <- ""
genes_celltypemarker_df$Cell_Type1[genes_celltypemarker_df$Cell_Type_Short %in% c("Distalconvolutedtubules")] <- "Distal convoluted tubule"
genes_celltypemarker_df$Cell_Type1[genes_celltypemarker_df$Cell_Type_Short %in% c("LoopofHenle")] <- "Loop of Henle"
genes_celltypemarker_df$Cell_Type1[genes_celltypemarker_df$Cell_Type_Short %in% c("Podocytes")] <- "Podocytes"
genes_celltypemarker_df$Cell_Type1[genes_celltypemarker_df$Cell_Type_Short %in% c("Proximaltubules")] <- "Proximal tubule"
genes_celltypemarker_df$Cell_Type1[genes_celltypemarker_df$Cell_Type_Short %in% c("Intercalatedcells", "PrincipleCells")] <- "Collecting duct"
genes_celltypemarker_df$Cell_Type1[genes_celltypemarker_df$Cell_Type_Short %in% c("Endothelial")] <- "Endothelial cells"
genes_celltypemarker_df$Cell_Type1[genes_celltypemarker_df$Cell_Type_Short %in% c("Fibroblasts", "Myofibroblasts")] <- "Fibroblasts"
genes_celltypemarker_df$Cell_Type1[genes_celltypemarker_df$Cell_Type_Short %in% c("Macrophages")] <- "Myleoid lineage immune cells"

## add cell type 2
genes_celltypemarker_df$Cell_Type2 <- ""
genes_celltypemarker_df$Cell_Type2[genes_celltypemarker_df$Cell_Type_Short %in% c("Intercalatedcells")] <- "Intercalated cells"
genes_celltypemarker_df$Cell_Type2[genes_celltypemarker_df$Cell_Type_Short %in% c("PrincipleCells")] <- "Principle cells"
genes_celltypemarker_df$Cell_Type2[genes_celltypemarker_df$Cell_Type_Short %in% c("Myofibroblasts")] <- "Myofibroblasts"
genes_celltypemarker_df$Cell_Type2[genes_celltypemarker_df$Cell_Type_Short %in% c("Macrophages")] <- "Monocytic lineage immune cells"

## add cell type 3
genes_celltypemarker_df$Cell_Type3 <- ""
genes_celltypemarker_df$Cell_Type3[genes_celltypemarker_df$Cell_Type_Short %in% c("Macrophages")] <- "Macrophages"
genes_celltypemarker_df$Cell_Type4 <- ""
## add species
genes_celltypemarker_df$Species <- "Mouse"

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "celltypemarkergenes_mouse.", run_id, ".tsv")
write.table(x = genes_celltypemarker_df, file = file2write, quote = F, sep = "\t", row.names = F)
