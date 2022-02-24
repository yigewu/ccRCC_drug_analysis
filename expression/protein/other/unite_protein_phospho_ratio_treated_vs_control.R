# Yige Wu @ WashU 2021 Jan

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
## input protein ratio treated vs. control
protein_ratio_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/get_ratio_treated_vs_control_bylog2protein/20220110.v1/log2Protein_ratio_treated_vs_control.20220110.v1.tsv")
## input phosphorylation ratio treated vs. control
phospho_ratio_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/phosphoprotein/get_foldchange/get_ratio_treated_vs_control_bylog2phospho/20220110.v1/Log2Phosphorylation_ratio_treated_vs_control.20220110.v1.tsv")

# combine -----------------------------------------------------------------
colnames_idx <- c("PG.Gene", "PG.ProteinName", "PG.ProteinDescriptions", "PG.ProteinAccession", "PTM_Name")
sample_ids <- colnames(phospho_ratio_df); sample_ids <- sample_ids[!(sample_ids %in% colnames_idx)]

protein_ratio_df$PTM_Name <- "Protein"
combinded_ratio_df <- rbind(protein_ratio_df[, c(colnames_idx, sample_ids)],
                            phospho_ratio_df[, c(colnames_idx, sample_ids)])

# calculate mean log2 ratio -----------------------------------------------
combinded_log2ratio_df <- cbind(combinded_ratio_df[,colnames_idx], log2(combinded_ratio_df[,sample_ids]))
combinded_log2ratio_df$mean_log2ratio_treatedvscontrol <- rowMeans(combinded_log2ratio_df[, sample_ids], na.rm = T)
combinded_log2ratio_df$number_nonna <- rowSums(!is.na(combinded_log2ratio_df[, sample_ids]))

combinded_log2ratio_filered_df <- combinded_log2ratio_df %>%
  filter(number_nonna >= )