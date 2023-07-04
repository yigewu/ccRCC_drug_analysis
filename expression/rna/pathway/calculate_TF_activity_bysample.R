# Yige Wu @ WashU 2023 Jun

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("OmnipathR")
packages = c(
  "plyr",
  "stringr",
  "reshape2",
  "data.table",
  "dplyr",
  "OmnipathR"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

# input dependencies ------------------------------------------------------
## input the gene expression TPM data
exp_df <- fread("./Data_Freeze/v1.dataFreeze.washU_rcc/3.geneExp/v3.20210116/datafreeze.v3.kallisto.geneExp.tsv", data.table = F)
## input which TFs to process
tf_target_df = import_transcriptional_interactions()
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/0.sample_info/v3.20210116/RCC_PDX_Samples.20210115.v2.xlsx")

# filter the TF-target by detection and direction -------------------------
nrow(tf_target_df)
tf_target_process_df = tf_target_df %>%
  filter(is_inhibition != 1) %>%
  filter(target_genesymbol %in% exp_df$Name)
rm(tf_target_df)
nrow(tf_target_process_df)
length(unique(tf_target_process_df$source_genesymbol)) ## 586 TFs

# preprocess expression ---------------------------------------------------
# summary(unlist(exp_df[,-1])) # TPM
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.0      0.0      0.1     17.0      1.8 376358.0 
## get selected gene expression sample ids
unique(sampleinfo_df$ShortTag)
sampleinfo_filtered_df = sampleinfo_df %>%
  filter(DataType == "RNA-Seq") %>%
  filter(ShortTag %in% c("Baseline", "Control", "Treated.Cabo", "Treated.Sap", "Treated.Cabo+Sap"))
log2tpm_df = log2(exp_df[,sampleinfo_filtered_df$Analysis_ID]+1)
rownames(log2tpm_df) = exp_df$Name
log2tpm_scale_df = t(apply(log2tpm_df, 1, scale))
colnames(log2tpm_scale_df) = colnames(log2tpm_df)

# process by TF -----------------------------------------------------------
tf_score_agg_df = NULL
for (tf_genesymbol in unique(tf_target_process_df$source_genesymbol)) {
  target_genesymbols = tf_target_process_df$target_genesymbol[tf_target_process_df$source_genesymbol == tf_genesymbol]
  # target_exp_df = log2tpm_df[target_genesymbols,]
  target_exp_df = log2tpm_scale_df[target_genesymbols,]
  if (length(target_genesymbols) > 1) {
    tf_score_vec = colMeans(target_exp_df)
  } else {
    tf_score_vec = target_exp_df
  }
  tf_score_df = data.frame(TF_genesymbol = tf_genesymbol, TF_activity_score = tf_score_vec, analysis_id = names(tf_score_vec))
  tf_score_agg_df = rbind(tf_score_agg_df, tf_score_df)
}
tf_score_wide_df = dcast(data = tf_score_agg_df, formula = TF_genesymbol ~ analysis_id, value.var = "TF_activity_score")

tf_score_agg_df = merge(x = tf_score_agg_df, y = sampleinfo_df, 
                        by.x = "analysis_id", by.y = "Analysis_ID", all.x = T)

# write -------------------------------------------------------------------
## set run id
version_tmp <- 3
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_drug_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

write.table(x = tf_score_agg_df, file = paste0(dir_out, "TF_activity.long.", run_id, ".tsv"), quote = F, row.names = F, sep = "\t")
write.table(x = tf_score_wide_df, file = paste0(dir_out, "TF_activity.wide.", run_id, ".tsv"), quote = F, row.names = F, sep = "\t")

