# Yige Wu @ WashU 2020 May

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
## input sample info to add in information about the model id, passage and TumorTissue info
b1_b8_samples_df1 <- readxl::read_excel(path = "./Resources/Bulk_Processed_Data/Sample_Info/b1_8/sampleInfo.washU_b1-b8.pdmr.other.passed.v5.extra.20200401.xlsx")
b1_b8_samples_df2 <- readxl::read_excel(path = "./Resources/Bulk_Processed_Data/Sample_Info/b1_8/sampleInfo.20200316.xlsx")

## inpu Xiaolu's PDX collection info to fill in the treatment timepont (1 month vs 2 month)
pdxcollection_df <- readxl::read_excel(path = "./Resources/RCC_model_info/RCC_PDX_Organoid_MasterLists/20200409_RCC_PDX_Organoid_MasterLists_XY.xlsx", sheet = "20200318_PDX sample collection")

#  merge and correct the sample info -------------------------------------------------------
## filter sample info 1
b1_b8_samples_df1 <- b1_b8_samples_df1 %>%
  filter(Source == "WashU") %>%
  filter(CancerGroup == "Kidney")
unique(b1_b8_samples_df1$Batch)
## get caseid to model id
caseid2modelid_df <- b1_b8_samples_df1 %>%
  filter(ModelID != "-") %>%
  select(Original_CaseID, ModelID) %>%
  unique()
## filter sample info 2
b1_b8_samples_df2 <- b1_b8_samples_df2 %>%
  filter(!grepl(x = Use, pattern = "No:"))
b1_b8_samples_df2$NCI_Passage[b1_b8_samples_df2$Group == "Human_Normal"] <- "HN"
b1_b8_samples_df2$NCI_Passage[b1_b8_samples_df2$Group == "Human_Tumor"] <- "HT"
b1_b8_samples_df2$SampleID[b1_b8_samples_df2$SampleID == "DRESL6C 3965M"] <- "DRESL6C_3965m"
b1_b8_samples_df2$SampleID[b1_b8_samples_df2$SampleID == "TRESL6C 3965M"] <- "TRESL6C_3965m"

## add model id to blood normal samples
b1_b8_samples_df2$ModelID[b1_b8_samples_df2$ModelID == "-"] <- mapvalues(x = b1_b8_samples_df2$Original_CaseID[b1_b8_samples_df2$ModelID == "-"], 
                                                                         from = caseid2modelid_df$Original_CaseID,
                                                                         to = caseid2modelid_df$ModelID)
b1_b8_samples_df1$ModelID[b1_b8_samples_df1$ModelID == "-"] <- mapvalues(x = b1_b8_samples_df1$Original_CaseID[b1_b8_samples_df1$ModelID == "-"], 
                                                                         from = caseid2modelid_df$Original_CaseID,
                                                                         to = caseid2modelid_df$ModelID)
## get common column names and merge
colnames_shared <- intersect(colnames(b1_b8_samples_df1), colnames(b1_b8_samples_df2))
colnames_shared
b1_b8_samples_df <- rbind(b1_b8_samples_df1[,colnames_shared], b1_b8_samples_df2[, colnames_shared])
## correct sample info
b1_b8_samples_df$SampleID[b1_b8_samples_df$Analysis_ID == "TWDE-WUR-062-DRESL6C_3965M"] <- "DRESL6C_3965m"
b1_b8_samples_df$TumorTissue[b1_b8_samples_df$TumorTissue == "PDX tumor;cabo treated"] <- "PDX tumor;treated"
## get unique
b1_b8_samples_df <- unique(b1_b8_samples_df)

# process treatment timepoint ---------------------------------------------
pdxcollection_df <- pdxcollection_df %>%
  mutate(SampleID.Xiaolu = gsub(x = `sample ID`, pattern = '[A-Z]', replacement = "")) %>%
  rename(TreatmentTimePeriod = ...6) %>%
  mutate(TreatmentTime.Start = str_split_fixed(string = TreatmentTimePeriod, pattern = "-", n = 2)[,1]) %>%
  mutate(TreatmentTime.Start = as.Date(TreatmentTime.Start, "%m/%d/%y")) %>%
  mutate(TreatmentTime.End = str_split_fixed(string = TreatmentTimePeriod, pattern = "-", n = 2)[,2]) %>%
  mutate(TreatmentTime.End = as.Date(TreatmentTime.End, "%m/%d/%y")) %>%
  mutate(Treatment.Days = TreatmentTime.End - TreatmentTime.Start)

# Batch 1- 8: summarize sequencing status ---------------------------------
b1_b8_pdx_df <- b1_b8_samples_df %>%
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
  mutate(WES = "Data Processed") %>%
  rename(Analysis_ID.WES = Analysis_ID) %>%
  select(col_names_keep, "WES", "Analysis_ID.WES")
b1_b8_rna_df <- b1_b8_pdx_df %>%
  filter(DataType == "RNA-Seq") %>%
  mutate(RNA = "Data Processed") %>%
  rename(Analysis_ID.RNA = Analysis_ID) %>%
  select(col_names_keep, "RNA", "Analysis_ID.RNA")
## add DNA sequecing status
b1_b8_seq_status_df <- merge(uniq_b1_b8_pdx_df, b1_b8_wes_df, by = col_names_keep, all.x = T)
## add RNA sequecing status
b1_b8_seq_status_df <- merge(b1_b8_seq_status_df, b1_b8_rna_df, by = col_names_keep, all.x = T)
## add the SampleID.Xiaolu and order
b1_b8_seq_status_df <- b1_b8_seq_status_df %>%
  mutate(SampleID.Xiaolu = ifelse(TumorTissue == "PDX tumor", gsub(x = str_split_fixed(string = SampleID.AcrossDataType, pattern = '_|\\-', n = 3)[,2], pattern = '[a-z]', replacement = "", ignore.case = T),
                                  ifelse(TumorTissue %in% c("PDX tumor;control", "PDX tumor;treated"), str_split_fixed(string = SampleID.AcrossDataType, pattern = "-", n = 3)[,2], NA))) %>%
  mutate(Treatment.Days = NA) %>%
  arrange(Formal_CaseID)
## add treatment days
b1_b8_seq_status_df$Treatment.Days = mapvalues(x = b1_b8_seq_status_df$SampleID.Xiaolu, from = pdxcollection_df$SampleID.Xiaolu, to = as.vector(pdxcollection_df$Treatment.Days))
b1_b8_seq_status_df$Treatment.Days[b1_b8_seq_status_df$Treatment.Days == b1_b8_seq_status_df$SampleID.Xiaolu] <- NA

# write output ------------------------------------------------------------
## write table
write.table(x = b1_b8_seq_status_df, file = paste0(dir_out, "RCC_PDX_Related_Samples.Data_Status.", "Batch1_8.", run_id, ".tsv"), quote = F, row.names = F, sep = "\t")


