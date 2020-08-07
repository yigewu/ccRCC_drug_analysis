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
batch1_8_data_status_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_status/write_batch1_8_data_status/20200612.v1/RCC_PDX_Related_Samples.Data_Status.Batch1_8.20200612.v1.tsv")
## input RCC PDX batch 9,
batch9_data_status_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_status/write_batch9_data_status/20200612.v1/RCC_PDX_Related_Samples.Data_Status.Batch9.20200612.v1.tsv")
## input RCC PDX batch 10
htan_b1_data_status_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_status/write_htan_batch1_data_status/20200731.v1/RCC_PDX_Related_Samples.Data_Status.HTAN_B1.20200731.v1.tsv")

# Merge data ---------------------------------
## merge batch1-8 with batch 10
model_seq_status_df <- merge(batch1_8_data_status_df, batch9_data_status_df, by = c("ModelID", "NCI_Passage", "SampleID.AcrossDataType", "Group", "TumorTissue",  "Batch", "WES", "RNA","SampleID.Xiaolu", "Treatment.Days", "Analysis_ID.WES", "Analysis_ID.RNA"), all = T)
## merge batch1-8 with batch 9
model_seq_status_df <- merge(model_seq_status_df, htan_b1_data_status_df, by = c("ModelID", "SampleID.AcrossDataType", "Group", "TumorTissue", "Batch", "WES", "RNA","SampleID.Xiaolu", "Treatment.Days", "Analysis_ID.WES", "Analysis_ID.RNA"), all = T)

# add treatment info ------------------------------------------------------
model_seq_status_df <- model_seq_status_df %>%
  mutate(Treated.Sun = grepl(x = SampleID.AcrossDataType, pattern = "Sun")) %>%
  mutate(Treated.Cab = grepl(x = SampleID.AcrossDataType, pattern = "Cab")) %>%
  mutate(Treated.Sap = grepl(x = SampleID.AcrossDataType, pattern = "Sap")) %>%
  mutate(Treated.ACF = grepl(x = SampleID.AcrossDataType, pattern = "ACF")) %>%
  mutate(Treated.Entinostat = grepl(x = SampleID.AcrossDataType, pattern = "Entinostat")) %>%
  mutate(TreatmentLengthGroup = ifelse(!is.na(Treatment.Days), 
                                       ifelse(abs(Treatment.Days - 30) < abs(Treatment.Days-60), "1_Month", "2_Months"), NA))

## order
model_seq_status_df <- model_seq_status_df %>%
  arrange(ModelID, Batch)

## motify the sequencing status
model_seq_status_df$WES[is.na(model_seq_status_df$WES)] <- "Not Sequenced"
model_seq_status_df$RNA[is.na(model_seq_status_df$RNA)] <- "Not Sequenced"

## write table
write.table(x = model_seq_status_df, file = paste0(dir_out, "RCC_PDX_Related_Samples.PDXNet_B1_9.HTAN_B1.Bulk_Data_Status.", run_id, ".tsv"), quote = F, row.names = F, sep = "\t")
