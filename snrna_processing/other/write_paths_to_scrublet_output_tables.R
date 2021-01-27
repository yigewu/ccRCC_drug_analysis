# Yige Wu @WashU Apr 2020
## filter the Cell Ranger filtered barcode by the same threshold
## create the seurat object

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
## input the cutoff summary
scrublet_cutoffs_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/Scrublet/summary/RCC_PDX_Scrublet_Cutoffs.csv")

# make table --------------------------------------------------------------
paths_df <- scrublet_cutoffs_df %>%
  filter(Optimal == "Yes") %>%
  mutate(Path_box = paste0(dir_base, "Resources/snRNA_Processed_Data/Scrublet/outputs/", Sample_id, "/Cutoff", Cutoff, "/", Sample_id, "_scrublet_output_table.csv.gz"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Paths_to_Scrulbet_Output_Tables.", run_id, ".tsv")
write.table(file = file2write, x = paths_df, quote = F, sep = "\t", row.names = F)

