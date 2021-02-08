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
## input the protein data
exp_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/preprocess/filter_and_transform_DIA_phosphorylation_data/20210205.v1/RCC_PDX.DIA_Phosphopeptide.Log2.Filtered.20210205.v1.tsv")
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")

# specify parameters --------------------------------------------------------------
group1_tmp <- "Con"
group2_tmp <- "Cabo"
group3_tmp <- "Sap"
group4_tmp <- "Cabo+ Sap"
min_datapoints <- 5
result_merged_df <- NULL
id_variables <- c("PG.ProteinAccessions", "PG.Genes", "PG.ProteinDescriptions", "PG.ProteinNames", "PG.ProteinGroups", "PEP.StrippedSequence")

# test each treatment group vs control ------------------------------------
## filter expression
exp_filtered_df <- exp_df[!duplicated(exp_df$PEP.StrippedSequence),]
exp_filtered_df <- exp_filtered_df %>%
  select(-EG.ModifiedSequence) %>%
  select(-EG.PrecursorId) %>%
  select(-PTM_Name_Aggr) %>%
  select(-PTM_Name_Human) %>%
  select(-PTM_Name_Mouse)
## get sample ids
ids_sample_group1 <- sampleinfo_df$`Sample ID`[sampleinfo_df$Treatment == group1_tmp]
ids_sample_group2 <- sampleinfo_df$`Sample ID`[sampleinfo_df$Treatment == group2_tmp]
ids_sample_group3 <- sampleinfo_df$`Sample ID`[sampleinfo_df$Treatment == group3_tmp]
ids_sample_group4 <- sampleinfo_df$`Sample ID`[sampleinfo_df$Treatment == group4_tmp]
## narrow down proteins to test
protein_group1_df <- exp_filtered_df[,ids_sample_group1]; datapoints_group1 <- rowSums(!is.na(protein_group1_df))
protein_group2_df <- exp_filtered_df[,ids_sample_group2]; datapoints_group2 <- rowSums(!is.na(protein_group2_df))
protein_group3_df <- exp_filtered_df[,ids_sample_group3]; datapoints_group3 <- rowSums(!is.na(protein_group3_df))
protein_group4_df <- exp_filtered_df[,ids_sample_group4]; datapoints_group4 <- rowSums(!is.na(protein_group4_df))
idxs_test <- which(datapoints_group1 >= min_datapoints & datapoints_group2 >= min_datapoints & datapoints_group3 >= min_datapoints & datapoints_group4 >= min_datapoints)
num_tests <- length(idxs_test)
## initiate
pvalue_vec <- vector(mode = "numeric", length = num_tests)
F_estimate_vec <- vector(mode = "numeric", length = num_tests)
## test per protein
for (i in 1:num_tests) {
  print(i)
  idx_test <- idxs_test[i]
  exp_tmp_df <- melt(data = exp_filtered_df[idx_test,], id.vars = id_variables)
  exp_tmp_df$Treatment <- mapvalues(x = exp_tmp_df$variable, from = sampleinfo_df$`Sample ID`, to = as.vector(sampleinfo_df$Treatment))
  # Compute the analysis of variance
  res.aov <- aov(value ~ Treatment, data = exp_tmp_df)
  # Summary of the analysis
  stat <- summary(res.aov)
  
  pvalue_vec[i] <- stat[[1]][["Pr(>F)"]][1]
  F_estimate_vec[i] <- stat[[1]][["F value"]][1]
}

## parse results
result_merged_df <- exp_df[idxs_test, id_variables]
result_merged_df <- cbind(result_merged_df,
                       data.frame(pvalue = pvalue_vec,
                                  F_estimate = F_estimate_vec))
result_merged_df$fdr <- p.adjust(p = result_merged_df$pvalue, method = "fdr")

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "Phosphosite.", "ANOVA.Treated_Groups_and_CT.", run_id, ".tsv")
write.table(x = result_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)
