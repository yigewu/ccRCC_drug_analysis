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
protein_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/preprocess/normalize_DIA_protein_data/20210111.v1/RCC_PDX.DIA_Protein.Log2.QuantileNormalized.20210111.v1.tsv")
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")

# specify parameters --------------------------------------------------------------
# unique(sampleinfo_df$Treatment)
group1_process <- c("Cabo", "Sap", "Cabo+ Sap")
group2_tmp <- "Con"
min_datapoints <- 5
result_merged_df <- NULL

# test each treatment group vs control ------------------------------------
# group1_tmp <- "Cabo"
for (group1_tmp in group1_process) {
  print(group1_tmp)
  ## get sample ids
  ids_sample_group1 <- sampleinfo_df$`Sample ID`[sampleinfo_df$Treatment == group1_tmp]
  ids_sample_group2 <- sampleinfo_df$`Sample ID`[sampleinfo_df$Treatment == group2_tmp]
  ## narrow down proteins to test
  protein_group1_df <- protein_df[,ids_sample_group1]
  datapoints_group1 <- rowSums(!is.na(protein_group1_df))
  protein_group2_df <- protein_df[,ids_sample_group2]
  datapoints_group2 <- rowSums(!is.na(protein_group2_df))
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
    protein_group1_tmp <- unlist(protein_group1_df[idx_test,])
    protein_group2_tmp <- unlist(protein_group2_df[idx_test,])
    ### test
    stat <- wilcox.test(x = protein_group1_tmp, y = protein_group2_tmp, conf.int = T)
    pvalue_vec[i] <- stat$p.value
    diff_estimate_vec[i] <- stat$estimate
    
    diff_bottom_vec[i] <- stat$conf.int[1]
    diff_top_vec[i] <- stat$conf.int[2]
    number_group1_vec[i] <- length(protein_group1_tmp[!is.na(protein_group1_tmp)])
    number_group2_vec[i] <- length(protein_group2_tmp[!is.na(protein_group2_tmp)])
  }
  
  ## parse results
  result_tmp_df <- protein_df[idxs_test, c("PG.ProteinAccessions", "PG.Genes", "PG.ProteinDescriptions", "PG.ProteinNames")]
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


# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "Wilcox.Each_Treated_Group_vs_CT.", run_id, ".tsv")
write.table(x = result_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)
