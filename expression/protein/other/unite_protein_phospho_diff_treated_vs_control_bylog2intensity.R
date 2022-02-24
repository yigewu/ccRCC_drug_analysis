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
protein_ratio_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/get_diff_treated_vs_control_bylog2protein/20220110.v1/log2Protein_diff_treated_vs_control.20220110.v1.tsv")
## input phosphorylation ratio treated vs. control
phospho_ratio_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/phosphoprotein/get_foldchange/get_diff_treated_vs_control_bylog2phospho/20220110.v1/Log2Phosphorylation_diff_treated_vs_control.20220110.v1.tsv")

# combine -----------------------------------------------------------------
colnames_idx <- c("PG.Gene", "PG.ProteinName", "PG.ProteinDescriptions", "PG.ProteinAccession", "PTM_Name")
sample_ids <- colnames(phospho_ratio_df); sample_ids <- sample_ids[!(sample_ids %in% colnames_idx)]

protein_ratio_df$PTM_Name <- "Protein"
combinded_ratio_df <- rbind(protein_ratio_df[, c(colnames_idx, sample_ids)],
                            phospho_ratio_df[, c(colnames_idx, sample_ids)])

# calculate mean log2 ratio -----------------------------------------------
sampleids_cabo <- sample_ids[grepl(pattern = "_Cabo_", x = sample_ids)]
combinded_ratio_df$mean_log2ratio_CabovsCon <- rowMeans(combinded_ratio_df[, sampleids_cabo], na.rm = T)
combinded_ratio_df$number_nonna_Cabo <- rowSums(!is.na(combinded_ratio_df[, sampleids_cabo]))
sampleids_sap <- sample_ids[grepl(pattern = "_Sap_", x = sample_ids)]
combinded_ratio_df$mean_log2ratio_SapvsCon <- rowMeans(combinded_ratio_df[, sampleids_sap], na.rm = T)
combinded_ratio_df$number_nonna_Sap <- rowSums(!is.na(combinded_ratio_df[, sampleids_sap]))
sampleids_cabosap <- sample_ids[grepl(pattern = "_Cabo\\+ Sap_", x = sample_ids)]
combinded_ratio_df$mean_log2ratio_CaboSapvsCon <- rowMeans(combinded_ratio_df[, sampleids_cabosap], na.rm = T)
combinded_ratio_df$number_nonna_CaaboSap <- rowSums(!is.na(combinded_ratio_df[, sampleids_cabosap]))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Protein_Phospho_Diff_Log2Intensity.Treated_vs_Control.", run_id, ".tsv")
write.table(file = file2write, x = combinded_ratio_df, quote = F, sep = "\t", row.names = F)
