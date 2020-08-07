# Yige Wu @ WashU 2020 March
## get the sequencing status for all the RCC related samples
## highlight what data type is missing (WES/RNA)

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies --------------------------------------------
## input bulk data status
## input all ccRCC samples data info
bulk_data_status_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_status/write_bulk_data_status/20200731.v1/RCC_PDX_Related_Samples.PDXNet_B1_9.HTAN_B1.Bulk_Data_Status.20200731.v1.tsv")
## input snRNA data status
sn_data_status_df <- readxl::read_xlsx(path = "./Resources/snRNA_Data_Generation/Sample_Collection/snRNA_Sample_Collection.xlsx")

# add sn data status ------------------------------------------------------
data_status_df <- bulk_data_status_df %>%
  mutate(snRNA = ifelse(SampleID.AcrossDataType %in% sn_data_status_df$SampleID.AcrossDataType, "Data Processed", "Not Sequenced"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "RCC_PDX_Related_Samples.PDXNet_B1_9.HTAN_B1.Bulk_SingleNuclei_Data_Status.", run_id, ".tsv")
write.table(x = data_status_df, file = file2write, quote = F, sep = "\t", row.names = F)
