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
## input RCC PDX batch 9 submission
batch9_samples_df <- readxl::read_excel(path = "./Resources/Bulk_Data_Generation/batch9/2020-01-27/20200127_Chen_Lab_DNA_RNA_Submission_Corrected.xlsx")
## input data info with analysis ids for batch 9 
datainfo_batch9_df <- readxl::read_excel(path = "./Resources/Bulk_Processed_Data/Data_Files/batch9/washU-batch9.dataInfo.202004.xlsx")
## input analysis ids for the RNA samples
ids_analysis_rna <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Data_Files/batch9/geneExp/human.pdx_transcriptome.sample", col.names = "Analysis_ID")
## inpu Xiaolu's PDX collection info to fill in the treatment timepont (1 month vs 2 month)
pdxcollection_df <- readxl::read_excel(path = "./Resources/RCC_model_info/RCC_PDX_Organoid_MasterLists/20200521_RCC_PDX_Organoid_MasterLists_XY.YW_edited.xlsx", sheet = "20200519_PDX sample collection")

# process treatment timepoint ---------------------------------------------
pdxcollection_df <- pdxcollection_df %>%
  mutate(SampleID.Xiaolu = gsub(x = `sample ID`, pattern = '[A-Z]', replacement = "")) %>%
  rename(TreatmentTimePeriod = `treatment duration`) %>%
  mutate(TreatmentTime.Start = str_split_fixed(string = TreatmentTimePeriod, pattern = "-", n = 2)[,1]) %>%
  mutate(TreatmentTime.Start = as.Date(TreatmentTime.Start, "%m/%d/%Y")) %>%
  mutate(TreatmentTime.End = str_split_fixed(string = TreatmentTimePeriod, pattern = "-", n = 2)[,2]) %>%
  mutate(TreatmentTime.End = as.Date(TreatmentTime.End, "%m/%d/%Y")) %>%
  mutate(Treatment.Days = TreatmentTime.End - TreatmentTime.Start)

# process the rna analysis ids --------------------------------------------
# ids_analysis_rna$Index_Sequence <- sapply(X = ids_analysis_rna$Analysis_ID, FUN = function(x) { str_split(string = x, pattern = '\\.')[[1]][length(str_split(string = x, pattern = '\\.')[[1]])]})
# ids_analysis_rna$SampleID <- sapply(X = ids_analysis_rna$Analysis_ID, FUN = function(x) { })

## just select for resl
ids_analysis_rna <- ids_analysis_rna %>%
  filter(grepl(pattern = "RESL", x = Analysis_ID, ignore.case = T)) %>%
  mutate(SampleID.Xiaolu = str_split_fixed(string = Analysis_ID, pattern = "-|_", n = 3)[,2]) %>%
  mutate(SampleID.Xiaolu = gsub(x = SampleID.Xiaolu, pattern = "[A-Z]|\\.", replacement = "", ignore.case = T))

# Batch 9----------------------------------------
## format submission info
batch9_pdx_df <- datainfo_batch9_df %>%
  filter(CancerType == "KIRC") %>%
  mutate(SampleID.AcrossDataType = gsub(x = SampleID, pattern = "DRESL|TRESL", replacement = "RESL")) %>%
  mutate(SampleID.AcrossDataType = gsub(x = SampleID.AcrossDataType, pattern = "-DNA|-RNA", replacement = "")) %>%
  mutate(ModelID = str_split_fixed(string = SampleID.AcrossDataType, pattern = "_|-", n = 2)[,1]) %>%
  mutate(ModelID = gsub(x = ModelID, pattern = '[A-Z]', replacement = "")) %>%
  mutate(ModelID = paste0("RESL", ModelID)) %>%
  mutate(SampleID.Xiaolu = str_split_fixed(string = SampleID.AcrossDataType, pattern = "-|_", n = 3)[,2]) %>%
  mutate(SampleID.Xiaolu = gsub(x = SampleID.Xiaolu, pattern = "[A-Z]|\\.", replacement = "", ignore.case = T))
## set column names to keep
batch9_col_names_keep <- c("ModelID", "SampleID.AcrossDataType", "NCI_Passage", "SampleID.Xiaolu")
## get unique info for each model
unique_batch9_pdx_df <- batch9_pdx_df %>%
  select(batch9_col_names_keep) %>%
  unique()
## filter sample by data type
batch9_wes_df <- batch9_pdx_df %>%
  filter(DataType == "WES") %>%
  mutate(WES = "Data Processed") %>%
  rename(Analysis_ID.WES = Analysis_ID) %>%
  select(batch9_col_names_keep, "WES", "Analysis_ID.WES")
## because we only have data info on WES, we are going to use the WES data info to fill in the RNA data info
batch9_rna_df <- batch9_pdx_df %>%
  filter(DataType == "WES") %>%
  mutate(DataType == "RNA-Seq") %>%
  mutate(RNA = "Data Processed") %>%
  select(batch9_col_names_keep, "RNA")
batch9_rna_df$Analysis_ID.RNA <- mapvalues(x = batch9_rna_df$SampleID.Xiaolu, from = ids_analysis_rna$SampleID.Xiaolu, to = ids_analysis_rna$Analysis_ID)
## add DNA sequecing status
batch9_seq_status_df <- merge(unique_batch9_pdx_df, batch9_wes_df, by = batch9_col_names_keep, all.x = T)
## add RNA sequecing status
batch9_seq_status_df <- merge(batch9_seq_status_df, batch9_rna_df, by = batch9_col_names_keep, all.x = T)
## add other info
batch9_seq_status_df <- batch9_seq_status_df %>%
  mutate(Group = "PDX") %>%
  mutate(Batch = "batch9") %>%
  mutate(TumorTissue = ifelse(grepl(x = SampleID.AcrossDataType, pattern = "-CT"), "PDX tumor;control", 
                              ifelse(grepl(x = SampleID.AcrossDataType, pattern = "Cab|Sap|Entinostat|ACF|Sun"), "PDX tumor;treated", "PDX tumor")))
## add treatment days
batch9_seq_status_df$Treatment.Days = mapvalues(x = batch9_seq_status_df$SampleID.Xiaolu, from = pdxcollection_df$SampleID.Xiaolu, to = as.vector(pdxcollection_df$Treatment.Days))
batch9_seq_status_df$Treatment.Days[batch9_seq_status_df$Treatment.Days == batch9_seq_status_df$SampleID.Xiaolu] <- NA

# write output ------------------------------------------------------------
write.table(x = batch9_seq_status_df, file = paste0(dir_out, "RCC_PDX_Related_Samples.Data_Status.", "Batch9.", run_id, ".tsv"), quote = F, row.names = F, sep = "\t")
