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
## input the barcode2celltype
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/snrna_processing/map_barcode/combine_human_and_mouse_barcode2celltype/20210225.v1/Human_and_Mouse.Barcode2CellType.20210225.v1.tsv", data.table = F)

# calculate fraction ------------------------------------------------------
barcodes_by_celltype_sample_df <- barcode2celltype_df %>%
  select(Cell_Type.Short, Id_Sample) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 0)
barcodes_by_sample_df <- barcode2celltype_df %>%
  select(Id_Sample) %>%
  table() %>%
  as.data.frame() %>%
  rename(Id_Sample = '.')
## merge number of barcodes
barcodes_by_celltype_sample_df <- merge(barcodes_by_celltype_sample_df, barcodes_by_sample_df, by = c("Id_Sample"), suffixes = c("_CellType_Sample", "_Sample"))
## calculate fraction
barcodes_by_celltype_sample_df <- barcodes_by_celltype_sample_df %>%
  mutate(Fraction_CellType_Sample = Freq_CellType_Sample/Freq_Sample)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "CellTypeFraction.", run_id, ".tsv")
write.table(x = barcodes_by_celltype_sample_df, file = file2write, sep = "\t", quote = F, row.names = F)

