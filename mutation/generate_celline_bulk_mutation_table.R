# Yige Wu @ WashU 2020 Feb

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

# input dependencies --------------------------------------------
## input mutation table from the data freeze
maf_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Data_Files/batch8_cellline/somaticMut/somaticMut_cellLine_hk2.dnp.meta2.non-silent.tsv")

# write mutation short amino acid change matrix ---------------------------------------
mut_mat <- get_somatic_mutation_detailed_matrix(pair_tab = ccRCC_SMGs, maf = maf_df)
write.table(x = mut_mat, file = paste0(dir_out, "RCC_CellLine.", "Mutation_Matrix.", run_id, ".tsv"), quote = F, row.names = F, sep = "\t")

# write mutation VAF matrix ---------------------------------------
vaf_mat <- get_somatic_mutation_vaf_matrix(pair_tab = ccRCC_SMGs, maf = maf_df)
write.table(x = vaf_mat, file = paste0(dir_out, "RCC_CellLine.", "Mutation_VAF_Matrix.", run_id, ".tsv"), quote = F, row.names = F, sep = "\t")

