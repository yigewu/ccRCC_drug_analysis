# Yige Wu @ WashU 2020 March
## get the sequencing status for all the RCC related samples
## highlight what data type is missing (WES/RNA/snRNA)

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
## input RCC PDX batch 10
batch10_samples_df <- readxl::read_excel(path = "./Resources/Bulk_Data_Generation/WES_RNA/HTAN_batch1//2020-03-04/HTAN_Batch1_Submission_Form_v2.xlsx")
## inpu Xiaolu's PDX collection info to fill in the treatment timepont (1 month vs 2 month)
pdxcollection_df <- readxl::read_excel(path = "./Resources/RCC_model_info/RCC_PDX_Organoid_MasterLists/20200521_RCC_PDX_Organoid_MasterLists_XY.YW_edited.xlsx", sheet = "20200519_PDX sample collection")

# process treatment timepoint ---------------------------------------------
pdxcollection_df <- pdxcollection_df %>%
  mutate(SampleID.Xiaolu = gsub(x = `sample ID`, pattern = '[A-Z]|#', replacement = "")) %>%
  mutate(SampleID.Xiaolu = str_split_fixed(string = SampleID.Xiaolu, pattern = " ", n = 2)[,1]) %>%
  rename(TreatmentTimePeriod = `treatment duration`) %>%
  mutate(TreatmentTime.Start.Raw = str_split_fixed(string = TreatmentTimePeriod, pattern = "-", n = 2)[,1]) %>%
  mutate(TreatmentTime.Start = as.Date(TreatmentTime.Start.Raw, "%m/%d/%Y")) %>%
  mutate(TreatmentTime.End.Raw = str_split_fixed(string = TreatmentTimePeriod, pattern = "-", n = 2)[,2]) %>%
  mutate(TreatmentTime.End = as.Date(TreatmentTime.End.Raw, "%m/%d/%Y")) %>%
  mutate(Treatment.Days = TreatmentTime.End - TreatmentTime.Start)

# Batch 10----------------------------------------
## filter to RESL
batch10_pdx_df <- batch10_samples_df %>%
  filter(`Tissue Name` == "PDX tumor") %>%
  rename(SampleName = `Sample Name`) %>%
  mutate(SampleID.AcrossDataType = str_split_fixed(string = SampleName, pattern = "-DNA|-RNA", n = 2)[,1]) %>%
  mutate(ModelID_tmp = str_split_fixed(string = `Individual Name`, pattern = "-", n = 2)[,1]) %>%
  mutate(ModelNo = gsub(x = ModelID_tmp, pattern = '[a-z]', replacement = "", ignore.case = T)) %>%
  mutate(ModelID = paste0("RESL", ModelNo))
## set column names to keep
batch10_col_names_keep <- c("ModelID", "SampleID.AcrossDataType")
## get unique info for each model
unique_batch10_pdx_df <- batch10_pdx_df %>%
  select(ModelID, SampleID.AcrossDataType) %>%
  unique()
## filter sample by data type
batch10_wes_df <- batch10_pdx_df %>%
  filter(`DNA Type` == "Genomic DNA") %>%
  mutate(WES = "Data Processed") %>%
  select(batch10_col_names_keep, "WES")
batch10_rna_df <- batch10_pdx_df %>%
  filter(`DNA Type` == "Total RNA") %>%
  mutate(RNA = "Sequencing") %>%
  select(batch10_col_names_keep, "RNA")
## add DNA sequecing status
batch10_seq_status_df <- merge(unique_batch10_pdx_df, batch10_wes_df, by = batch10_col_names_keep, all.x = T)
## add RNA sequecing status
batch10_seq_status_df <- merge(batch10_seq_status_df, batch10_rna_df, by = batch10_col_names_keep, all.x = T)
## add other info
batch10_seq_status_df <- batch10_seq_status_df %>%
  mutate(Group = "PDX") %>%
  mutate(Batch = "HTAN_batch1") %>%
  mutate(TumorTissue = ifelse(grepl(x = SampleID.AcrossDataType, pattern = "-CT"), "PDX tumor;control", "PDX tumor;treated")) %>%
  mutate(SampleID.Xiaolu = str_split_fixed(string = SampleID.AcrossDataType, pattern = "-", n = 3)[,2]) %>%
  mutate(Analysis_ID.WES = NA) %>%
  mutate(Analysis_ID.RNA = NA)

## add treatment days
batch10_seq_status_df$Treatment.Days = mapvalues(x = batch10_seq_status_df$SampleID.Xiaolu, from = pdxcollection_df$SampleID.Xiaolu, to = as.vector(pdxcollection_df$Treatment.Days))

# write output ------------------------------------------------------------
write.table(x = batch10_seq_status_df, file = paste0(dir_out, "RCC_PDX_Related_Samples.Data_Status.", "HTAN_B1.", run_id, ".tsv"), quote = F, row.names = F, sep = "\t")

