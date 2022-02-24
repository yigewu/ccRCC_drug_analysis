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
## input the protein z-score
protein_zscore_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/get_protein_zscore_acrossamples/20220110.v1/Protein_zscore.20220110.v1.tsv")
## input the phosphorylation z-score
phospho_zscore_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/phosphoprotein/other/get_phosphorylation_zscore_acrossamples/20220110.v1/Phosphorylation_zscore.20220110.v1.tsv")
## input the sample info
meta_data_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")

# calculate by pathway ----------------------------------------------------
sample_ids <- meta_data_df$`Sample ID`
pathwaynames_process <- c("Apoptosis", "PI3K_Akt", "TSC_mTOR")
pathway_tmp <- pathwaynames_process[3]
pathwayscores_df <- NULL
pathwaynames_vec <- NULL
for (pathway_tmp in pathwaynames_process) {
  ## input the pathway protein list
  pathway2members_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Gene_Lists/Pathway_score_members.xlsx", sheet = pathway_tmp)
  pathway2protein_df <- pathway2members_df  %>%
    filter(Is_phospho == "No")
  pathway2phospho_df <- pathway2members_df  %>%
    filter(Is_phospho == "Yes") %>%
    as.data.frame()
  if (nrow(pathway2protein_df) > 0) {
    exp_process_df <- merge(x = protein_zscore_df, y = pathway2protein_df, by.x = c("PG.Gene"), by.y = c("Gene_symbol"))
    exp_positive_df <- exp_process_df %>%
      filter(Type_of_regulation == "Positive")
    if (nrow(exp_negative_df) > 0) {
      pathwayscores_pos_vec <- colMeans(exp_positive_df[,sample_ids], na.rm = T)
    } else {
      pathwayscores_pos_vec <- 0
    }
    exp_negative_df <- exp_process_df %>%
      filter(Type_of_regulation == "Negative")
    if (nrow(exp_negative_df) > 0) {
      pathwayscores_neg_vec <- colMeans(exp_negative_df[,sample_ids], na.rm = T)
    } else {
      pathwayscores_neg_vec <- 0
    }
    pathwayscores_protein_vec <- (pathwayscores_pos_vec - pathwayscores_neg_vec)
  } else {
    pathwayscores_protein_vec <- rep(0, length(sample_ids))
  }
  if (nrow(pathway2phospho_df) > 0) {
    exp_process_df <- merge(x = phospho_zscore_df, y = pathway2phospho_df, by.x = c("PG.Gene", "PTM_Name"), by.y = c("Gene_symbol", "Site_phospho"))
    if (nrow(exp_process_df) > 0) {
      exp_positive_df <- exp_process_df %>%
        filter(Type_of_regulation == "Positive")
      if (nrow(exp_negative_df) > 0) {
        pathwayscores_pos_vec <- colMeans(exp_positive_df[,sample_ids], na.rm = T)
      } else {
        pathwayscores_pos_vec <- 0
      }
      exp_negative_df <- exp_process_df %>%
        filter(Type_of_regulation == "Negative")
      if (nrow(exp_negative_df) > 0) {
        pathwayscores_neg_vec <- colMeans(exp_negative_df[,sample_ids], na.rm = T)
      } else {
        pathwayscores_neg_vec <- 0
      }
      pathwayscores_phospho_vec <- (pathwayscores_pos_vec - pathwayscores_neg_vec)
    } else {
      pathwayscores_phospho_vec <- rep(0, length(sample_ids))
    }
  }
  if (all(pathwayscores_phospho_vec == 0) & all(pathwayscores_protein_vec == 0)) {
    cat(paste0(pathway_tmp, ": No pathway protein/phosphosite detected!\n"))
  } else {
    pathwayscores_vec <- colMeans(rbind(pathwayscores_protein_vec, pathwayscores_phospho_vec), na.rm = T)
    pathwayscores_df <- rbind(pathwayscores_df, matrix(pathwayscores_vec, nrow = 1))
    pathwaynames_vec <- c(pathwaynames_vec, pathway_tmp)
  }
}
pathwayscores_df <- data.frame(pathwayscores_df)
colnames(pathwayscores_df) <- sample_ids
pathwayscores_df <- cbind(data.frame(Pathway_name = pathwaynames_vec), pathwayscores_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Pathway_scores_bysample.", run_id, ".tsv")
write.table(x = pathwayscores_df, file = file2write, sep = "\t", row.names = F, quote = F)
