# Yige Wu @ WashU 2021 Jan

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
packages_defaultinstall <- c("glmnet", "itertools", "missForest", "randomForest", "cluster", "survival", "Rcpp", "foreach","iterators", "Matrix", "devtools")
for (pkg_name_tmp in packages_defaultinstall) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  library(package = pkg_name_tmp, character.only = T)
}
if (!("impute" %in% installed.packages()[,1])) {
  BiocManager::install("impute")
}
library(package = "impute", character.only = T)
require("remotes")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS= "true")
install_github("WangLab-MSSM/DreamAI/Code")
library(DreamAI)
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/load_pkgs.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the protein data
# protein_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/preprocess/normalize_DIA_protein_data/20210111.v1/RCC_PDX.DIA_Protein.Log2.QuantileNormalized.20210111.v1.tsv")
protein_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/preprocess/filter_and_transform_DIA_protein_data/20210205.v1/RCC_PDX.DIA_Protein.Log2.Filtered.20210205.v1.tsv")

# run impute --------------------------------------------------------------
## ADMIN, RegImpute, MissForest, KNN work
## SpectroFM throws error, Birnn does not work
## SpectroFM: Error in do_libfm(strvectors, args) : stof: no conversion
regimpute_result<- DreamAI(data = protein_df[, 5:ncol(protein_df)],
                 k=10,maxiter_MF = 10, ntree = 100,maxnodes = NULL,maxiter_ADMIN=30,tol=10^(-2),gamma_ADMIN=NA,gamma=50,CV=FALSE,fillmethod="row_mean",maxiter_RegImpute=10,conv_nrmse = 1e-6,iter_SpectroFM=40, 
                 method = c("RegImpute"),out="Ensemble")
protein_regimputed_value_df <- regimpute_result$Ensemble
View(protein_regimputed_value_df)
protein_regimputed_df <- cbind(protein_df[, 1:4], protein_regimputed_value_df)

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "RCC_PDX.DIA_Protein.Log2.", "RegImpute", ".tsv")
write.table(x = protein_regimputed_df, file = file2write, quote = F, sep = "\t", row.names = F)
