# Yige Wu @ WashU 2020 Jul

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
library(readxl)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

#  input dependencies -----------------------------------------------------
## input bulk & sn data status
data_status_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_status/write_sn_bulk_data_status/20200731.v1/RCC_PDX_Related_Samples.PDXNet_B1_9.HTAN_B1.Bulk_SingleNuclei_Data_Status.20200731.v1.tsv")
## input the JHU samples
jhu_samples_df <- read_excel(path = "./Resources/RCC_model_info/PDX Sample Collection Tracking_20200731.xlsx", sheet = "To JHU ")

# get the number id for the jhu samples -----------------------------------
jhu_samples_df <- jhu_samples_df %>%
  filter(`Left tumor weight (mg)` != "L-piece 1" | is.na(`Left tumor weight (mg)`)) %>%
  mutate(SampleID.Xiaolu = gsub(x = str_split_fixed(string = `PDX line_Batch_ID`, pattern = "_", n = 3)[,3], pattern = "[A-Z]", replacement = ""))

# merge -------------------------------------------------------------------
data_status_filtered_df <- data_status_df %>%
  filter(!is.na(SampleID.Xiaolu))
jhu_samples_df$WES <- mapvalues(x = jhu_samples_df$SampleID.Xiaolu, from = data_status_filtered_df$SampleID.Xiaolu, to = as.vector(data_status_filtered_df$WES))
jhu_samples_df$WES[jhu_samples_df$WES == jhu_samples_df$SampleID.Xiaolu & jhu_samples_df$WES != ""] <- "Not Sequenced"
jhu_samples_df$RNA_Seq <- mapvalues(x = jhu_samples_df$SampleID.Xiaolu, from = data_status_filtered_df$SampleID.Xiaolu, to = as.vector(data_status_filtered_df$RNA))
jhu_samples_df$RNA_Seq[jhu_samples_df$RNA_Seq == jhu_samples_df$SampleID.Xiaolu & jhu_samples_df$RNA_Seq != ""] <- "Not Sequenced"
jhu_samples_df$snRNA_Seq <- mapvalues(x = jhu_samples_df$SampleID.Xiaolu, from = data_status_filtered_df$SampleID.Xiaolu, to = as.vector(data_status_filtered_df$snRNA))
jhu_samples_df$snRNA_Seq[jhu_samples_df$snRNA_Seq == jhu_samples_df$SampleID.Xiaolu & jhu_samples_df$snRNA_Seq != ""] <- "Not Sequenced"

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "PDX_Sample_Collection_Tracking.", "Sequence_Status_Annotated.", run_id, ".tsv")
write.table(x = jhu_samples_df, file = file2write, quote = F, sep = "\t", row.names = F)

