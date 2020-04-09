# Yige Wu @ WashU 2020 March
## get the sequencing status for all the RCC related samples
## highlight what data type is missing (WES/RNA/snRNA)

# set up libraries and output directory -----------------------------------
## set up working directory and source functions and load libraries
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/ccRCC_Drug_shared.R")
## set run id 
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies --------------------------------------------
## input sample info to add in information about the model id, passage and TumorTissue info
b1_b8_samples_df <- readxl::read_excel(path = "./PDX-Pilot/DataFreeze/Sample_info/sampleInfo.v4/sampleInfo.washU_b1-b8.pdmr.other.passed.v4.20200302.xlsx")
## input RCC PDX batch 9
batch9_samples_df <- readxl::read_excel(path = "./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/Resources/Bulk_Data_Generation/2020-01-27/20200127_Chen_Lab_DNA_RNA_Submission_Corrected.xlsx")
## input RCC PDX batch 10
batch10_samples_df <- readxl::read_excel(path = "./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/Resources/Bulk_Data_Generation/2020-03-04/HTAN_Batch1_Submission_Form_v2.xlsx")

# Batch 1- 8: summarize sequencing status ---------------------------------
## filter to kidney
b1_b8_pdx_df <- b1_b8_samples_df %>%
  filter(Source == "WashU") %>%
  filter(CancerGroup == "Kidney") %>%
  mutate(SampleID.AcrossDataType = gsub(x = SampleID, pattern = "DRESL|TRESL", replacement = "RESL")) %>%
  mutate(SampleID.AcrossDataType = gsub(x = SampleID.AcrossDataType, pattern = "-DNA|-RNA", replacement = ""))
b1_b8_pdx_df$ModelID[b1_b8_pdx_df$ModelID == "-"] <- plyr::mapvalues(x = b1_b8_pdx_df$Formal_CaseID[b1_b8_pdx_df$ModelID == "-"], from = b1_b8_pdx_df$Formal_CaseID[b1_b8_pdx_df$Group == "PDX"], to = as.vector(b1_b8_pdx_df$ModelID[b1_b8_pdx_df$Group == "PDX"]))
## set column names to keep
col_names_keep <- c("Original_CaseID", "Formal_CaseID",  "Group", "ModelID", "TumorTissue", "NCI_Passage", "Batch", "SampleID.AcrossDataType")
## get unique info for each model
uniq_b1_b8_pdx_df <- b1_b8_pdx_df %>%
  select(col_names_keep) %>%
  unique() %>%
  arrange(Formal_CaseID)
## filter sample by data type
b1_b8_wes_df <- b1_b8_pdx_df %>%
  filter(DataType == "WES") %>%
  mutate(WES = "Data Received") %>%
  select(col_names_keep, "WES")
b1_b8_rna_df <- b1_b8_pdx_df %>%
  filter(DataType == "RNA-Seq") %>%
  mutate(RNA = "Data Received") %>%
  select(col_names_keep, "RNA")
## add DNA sequecing status
b1_b8_seq_status_df <- merge(uniq_b1_b8_pdx_df, b1_b8_wes_df, by = col_names_keep, all.x = T)
## add RNA sequecing status
b1_b8_seq_status_df <- merge(b1_b8_seq_status_df, b1_b8_rna_df, by = col_names_keep, all.x = T)
## order
b1_b8_seq_status_df <- b1_b8_seq_status_df %>%
  arrange(Formal_CaseID)

# Batch 9----------------------------------------
## filter to RESL
batch9_pdx_df <- batch9_samples_df %>%
  filter(grepl(x = `Sample ID`, pattern = "RESL")) %>%
  rename(SampleID = `Sample ID`) %>%
  rename(ModelID = `PDX ID`) %>%
  mutate(SampleID.AcrossDataType = str_split_fixed(string = SampleID, pattern = "-DNA|-RNA", n = 2)[,1])
## set column names to keep
batch9_col_names_keep <- c("ModelID", "SampleID.AcrossDataType")
## get unique info for each model
unique_batch9_pdx_df <- batch9_pdx_df %>%
  select(ModelID, SampleID.AcrossDataType) %>%
  unique()
## filter sample by data type
batch9_wes_df <- batch9_pdx_df %>%
  filter(`Sample Type` == "genomic DNA") %>%
  mutate(WES = "Sequencing") %>%
  select(batch9_col_names_keep, "WES")
batch9_rna_df <- batch9_pdx_df %>%
  filter(`Sample Type` == "total RNA") %>%
  mutate(RNA = "Sequencing") %>%
  select(batch9_col_names_keep, "RNA")
## add DNA sequecing status
batch9_seq_status_df <- merge(unique_batch9_pdx_df, batch9_wes_df, by = batch9_col_names_keep, all.x = T)
## add RNA sequecing status
batch9_seq_status_df <- merge(batch9_seq_status_df, batch9_rna_df, by = batch9_col_names_keep, all.x = T)
## add other info
batch9_seq_status_df <- batch9_seq_status_df %>%
  mutate(Group = "PDX") %>%
  mutate(Batch = "batch9") %>%
  mutate(TumorTissue = ifelse(grepl(x = SampleID.AcrossDataType, pattern = "-CT"), "PDX tumor;control", "PDX tumor;treated"))

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
  mutate(WES = "Sequencing") %>%
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
  mutate(Batch = "batch10") %>%
  mutate(TumorTissue = ifelse(grepl(x = SampleID.AcrossDataType, pattern = "-CT"), "PDX tumor;control", "PDX tumor;treated"))

# Merge data ---------------------------------
## merge batch1-8 with batch 10
model_seq_status_df <- merge(b1_b8_seq_status_df, batch9_seq_status_df, by = c("ModelID", "SampleID.AcrossDataType", "Group", "TumorTissue", "Batch", "WES", "RNA"), all = T)
## merge batch1-8 with batch 9
model_seq_status_df <- merge(model_seq_status_df, batch10_seq_status_df, by = c("ModelID", "SampleID.AcrossDataType", "Group", "TumorTissue", "Batch", "WES", "RNA"), all = T)

## order
model_seq_status_df <- model_seq_status_df %>%
  select(col_names_keep, "WES", "RNA") %>%
  arrange(ModelID)
## write table
write.table(x = model_seq_status_df, file = paste0(dir_out, "RCC_PDX_Related_Samples.Sequencing_Status.", run_id, ".tsv"), quote = F, row.names = F, sep = "\t")
