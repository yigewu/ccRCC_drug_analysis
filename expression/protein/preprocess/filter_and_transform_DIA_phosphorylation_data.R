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
ptm_values_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Data_Files/Protein/01252021/LiDing_PDX_ccRCC_204phospho_library_Report.txt")
## input phosphosite location
ptm_names_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/preprocess/aggregate_phosphosite_names/20210204.v1/Phosphorylation_Sites.Aggregated.20210204.v1.tsv")
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Resources/Bulk_Processed_Data/Data_Files/Protein/01062021/WUSTL to JHU_ccRCC PDX_sample information.xlsx")

# filter rows by phospho-peptide and get column names------------------------------------------
phospho_df <- ptm_values_df %>%
  filter(grepl(x = EG.ModifiedSequence, pattern = "Phospho"))
colnames_all <- colnames(phospho_df)
colnames_id <- colnames_all[1:8]
colnames_data <- colnames_all[grepl(pattern = "MS2Quantity", x = colnames_all)]

# map phosphosite names, change data column names and log2 -----------------------------------------------------
exp_data_df <- phospho_df[, colnames_data]
exp_id_df <- phospho_df[, colnames_id]

## map phosphosite names
exp_id_df <- merge(x = exp_id_df, y = ptm_names_df, by = c("PG.ProteinGroups", "PG.ProteinNames", "EG.ModifiedSequence"), all.x = T)

## create column name mapping result
colnames_data_df <- data.frame(column_name = colnames_data)
colnames_data_df <- colnames_data_df %>%
  mutate(sample_position = str_split_fixed(string = column_name, pattern = "JHU_LM2_LiDing_PDX_", n = 2)[,2]) %>%
  mutate(sample_position = str_split_fixed(string = sample_position, pattern = "\\.raw", n = 2)[,1])
## add sample id
colnames_data_df$sample_id <- plyr::mapvalues(x = colnames_data_df$sample_position, from = sampleinfo_df$`my position`, to = as.vector(sampleinfo_df$`Sample ID`))
## change column name
colnames_data_new <- as.vector(colnames_data_df$sample_id)
colnames(exp_data_df) <- colnames_data_new
## log2 data
exp_data_log2_df <- log2(x = exp_data_df)
## combine data
exp_log2_df <- cbind(exp_id_df, exp_data_log2_df)

# filter by data points ---------------------------------------------------
## filter by number of data points
counts_datapoints <- rowSums(!is.na(exp_log2_df[, colnames_data_new]))
exp_log2_filtered_df <- exp_log2_df[counts_datapoints >= 0.7*length(colnames_data_new),]

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "RCC_PDX.DIA_Phosphopeptide.Log2.", run_id, ".tsv")
write.table(x = exp_log2_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "RCC_PDX.DIA_Phosphopeptide.Log2.Filtered.", run_id, ".tsv")
write.table(x = exp_log2_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)

