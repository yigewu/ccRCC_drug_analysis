# Yige Wu @ WashU 2022 Jan

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
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
## input protein ratio treated vs. control
protein_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/preprocess/split_protein_data_bygene/20220120.v1/RCC_PDX.DIA_Protein.Log2.Bygene.20220120.v1.tsv")

## input phosphorylation ratio treated vs. control
phospho_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/phosphoprotein/preprocess/split_phosphorylation_bygene/20220120.v1/RCC_PDX.DIA_Phosphopeptide.Log2.Bygene.20220120.v1.tsv")

# combine -----------------------------------------------------------------
colnames_idx <- c("PG.Gene", "PG.ProteinName", "PG.ProteinDescriptions", "PG.ProteinAccession", "PTM_Name")
sample_ids <- colnames(phospho_df); sample_ids <- sample_ids[grepl(pattern = "RESL", x = sample_ids)]; sample_ids

protein_df$PTM_Name <- "Protein"
combinded_ratio_df <- rbind(protein_df[, c(colnames_idx, sample_ids)],
                            phospho_df[, c(colnames_idx, sample_ids)])
combinded_ratio_df$ID <- paste0(combinded_ratio_df$PG.Gene, "_", combinded_ratio_df$PTM_Name)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Protein_Phospho_Log2Intensity.", run_id, ".tsv")
write.table(file = file2write, x = combinded_ratio_df, quote = F, sep = "\t", row.names = F)
