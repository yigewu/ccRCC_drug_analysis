# Yige Wu @ WashU 2022 Feb

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
packages = c(
  "plyr",
  "stringr",
  "reshape2",
  "data.table",
  "dplyr"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

# input dependencies ------------------------------------------------------
## input the protein data
exp_df <- fread("./Data_Freeze/v1.dataFreeze.washU_rcc/3.geneExp/v3.20210116/datafreeze.v3.kallisto.geneExp.tsv", data.table = F)
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/0.sample_info/v3.20210116/RCC_PDX_Samples.20210115.v2.xlsx")

# specify parameters --------------------------------------------------------------
# unique(sampleinfo_df$Treatment)
group1_process <- c("Treated.Cabo+Sap")
group2_process <- c("Treated.Cabo", "Treated.Sap", "Control")
min_datapoints <- 5
result_merged_df <- NULL

# test each treatment group vs control ------------------------------------
## pre-process
sampleinfo_df <- sampleinfo_df %>%
  filter(DataType == "RNA-Seq") %>%
  mutate(Treatment = ShortTag) %>%
  mutate(Treatment_length = paste0(Treatment.Month, " month")) %>%
  mutate(PairTag.middle = str_split_fixed(string = Analysis_ID, pattern = "_", n = 3)[,2]) %>%
  mutate(PairTag = ifelse(!is.na(PairTag), PairTag, paste0(ModelID, "_", PairTag.middle))) %>%
  mutate(id_model_length = paste0(PairTag, "_", Treatment_length))
colnames(exp_df)

for (group1_tmp in group1_process) {
  print(group1_tmp)
  ## get sample ids
  sampleinfo_group1_df <- sampleinfo_df %>%
    filter(Treatment_length == "1 month") %>%
    filter(Treatment == group1_tmp)
  ids_model_length <- sampleinfo_group1_df$id_model_length
  ids_sample_group1 <- sampleinfo_group1_df$Analysis_ID
  
  ## narrow down proteins to test
  exp_group1_df <- exp_df[,ids_sample_group1]
  datapoints_group1 <- rowSums(!is.na(exp_group1_df))
  
  for (group2_tmp in group2_process) {
    sampleinfo_group2_df <- sampleinfo_df %>%
      filter(Treatment_length == "1 month") %>%
      filter(Treatment == group2_tmp)
    ids_sample_group2 <- mapvalues(x = ids_model_length, from = sampleinfo_group2_df$id_model_length, to = as.vector(sampleinfo_group2_df$Analysis_ID))
    
    exp_group2_df <- exp_df[,ids_sample_group2]
    datapoints_group2 <- rowSums(!is.na(exp_group2_df))
    idxs_test <- which(datapoints_group1 >= min_datapoints & datapoints_group2 >= min_datapoints)
    num_tests <- length(idxs_test)
    ## initiate
    pvalue_vec <- vector(mode = "numeric", length = num_tests)
    diff_estimate_vec <- vector(mode = "numeric", length = num_tests)
    diff_bottom_vec <- vector(mode = "numeric", length = num_tests)
    diff_top_vec <- vector(mode = "numeric", length = num_tests)
    number_group1_vec <- vector(mode = "numeric", length = num_tests)
    number_group2_vec <- vector(mode = "numeric", length = num_tests)
    ## test per protein
    for (i in 1:num_tests) {
      print(i)
      idx_test <- idxs_test[i]
      exp_group1_tmp <- unlist(exp_group1_df[idx_test,])
      exp_group2_tmp <- unlist(exp_group2_df[idx_test,])
      ### test
      stat <- t.test(x = exp_group1_tmp, y = exp_group2_tmp, conf.int = T, paired = T)
      pvalue_vec[i] <- stat$p.value
      diff_estimate_vec[i] <- stat$estimate
      
      diff_bottom_vec[i] <- stat$conf.int[1]
      diff_top_vec[i] <- stat$conf.int[2]
      number_group1_vec[i] <- length(exp_group1_tmp[!is.na(exp_group1_tmp)])
      number_group2_vec[i] <- length(exp_group2_tmp[!is.na(exp_group2_tmp)])
    }
    
    ## parse results
    result_tmp_df <- data.frame(Name = exp_df[idxs_test, c("Name")])
    result_tmp_df <- cbind(result_tmp_df,
                           data.frame(pvalue = pvalue_vec,
                                      diff_estimate = diff_estimate_vec,
                                      diff_bottom = diff_bottom_vec,
                                      diff_top = diff_top_vec,
                                      number_group1 = number_group1_vec,
                                      number_group2 = number_group2_vec,
                                      group1 = group1_tmp,
                                      group2 = group2_tmp))
    result_tmp_df$fdr <- p.adjust(p = result_tmp_df$pvalue, method = "fdr")
    result_merged_df <- rbind(result_merged_df, result_tmp_df)
  }
}

# save output -------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_drug_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
file2write <- paste0(dir_out, "mRNA.Ttest.Paired.1month.Combo_vs_single.", run_id, ".tsv")
write.table(x = result_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)
