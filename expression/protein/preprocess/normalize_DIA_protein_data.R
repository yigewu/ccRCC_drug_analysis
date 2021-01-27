# Yige Wu @ WashU 2021 Jan

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
if (!("preprocessCore" %in% installed.packages()[,1])) {
  #install if necessary
  BiocManager::install("preprocessCore")
}
#load package
library(preprocessCore)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the protein data
protein_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Data_Files/Protein/01062021/LM_LiDing_PDX_Report.xls", na.strings = "Filtered")
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Resources/Bulk_Processed_Data/Data_Files/Protein/01062021/WUSTL to JHU_ccRCC PDX_sample information.xlsx")

# change column names -----------------------------------------------------
colnames_all <- colnames(protein_df)
colnames_id <- colnames_all[1:4]
colnames_data <- colnames_all[5:length(colnames_all)]
## create column name mapping result
colnames_data_df <- data.frame(column_name = colnames_data)
colnames_data_df <- colnames_data_df %>%
  mutate(sample_position = str_split_fixed(string = column_name, pattern = "JHU_LM2_LiDing_PDX_", n = 2)[,2]) %>%
  mutate(sample_position = gsub(x = sample_position, pattern = '\\.raw\\.PG\\.Quantity', replacement = ""))
## add sample id
colnames_data_df$sample_id <- plyr::mapvalues(x = colnames_data_df$sample_position, from = sampleinfo_df$`my position`, to = as.vector(sampleinfo_df$`Sample ID`))
## change column name
colnames_data_new <- as.vector(colnames_data_df$sample_id)
protein_new_df <- protein_df
colnames(protein_new_df) <- c(colnames_id, colnames_data_new)

# log2 data ---------------------------------------------------------------
protein_data_df <- protein_new_df[, colnames_data_new]
protein_data_log2_df <- log2(x = protein_data_df)

# quantile normalization --------------------------------------------------
protein_data_log2_qn_mat <- preprocessCore::normalize.quantiles(x = as.matrix(protein_data_log2_df))
protein_data_log2_qn_df <- as.data.frame(protein_data_log2_qn_mat)
colnames(protein_data_log2_qn_df) <- colnames_data_new

# combine data ------------------------------------------------------------
protein_log2_qn_df <- cbind(protein_new_df[,colnames_id], protein_data_log2_qn_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "RCC_PDX.DIA_Protein.Log2.QuantileNormalized.", run_id, ".tsv")
write.table(x = protein_log2_qn_df, file = file2write, quote = F, sep = "\t", row.names = F)


