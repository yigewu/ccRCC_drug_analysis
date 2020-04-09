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
## input formated mouse marker gene list
gene2celltype_mouse_df <- fread(input = "./Resources/Analysis_Results/dependencies/format_ruiyang_mouse_marker_genes/20200409.v1/celltypemarkergenes_mouse.20200409.v1.tsv", data.table = F)
## input human marker gene list
gene2celltype_human_df <- fread(input = "./Resources/Gene_Lists/Cell_Type_Marker_Genes/RCC_Marker_Gene_Table - Gene2CellType_Tab.20200406.v1.tsv", data.table = F)

# rename columns ----------------------------------------------------------
gene2celltype_human_df <- gene2celltype_human_df %>%
  rename(Gene_Symbol = Gene) %>%
  mutate(Species = "Human")

# merge -------------------------------------------------------------------
colnames_shared <- intersect(colnames(gene2celltype_human_df), colnames(gene2celltype_mouse_df))
colnames_shared
gene2celltype_df <- rbind(gene2celltype_mouse_df[, colnames_shared],
                          gene2celltype_human_df[, colnames_shared])

# write output -------------------------------------------------------------
file2write <- paste0(dir_out, "celltypemarkergenes_mouse_human.rcc.", run_id, ".tsv")
write.table(x = gene2celltype_df, file = file2write, quote = F, sep = "\t", row.names = F)
