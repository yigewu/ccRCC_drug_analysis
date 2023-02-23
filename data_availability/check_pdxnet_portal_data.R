# Yige Wu @ WashU 2021 Nov

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
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

# input dependencies ------------------------------------------------------
## input sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/0.sample_info/v3.20210116/RCC_PDX_Samples.20210115.v2.xlsx")
## input the list of sample currently on PDXNet portal
dnarna_samples_df <- read_excel(path = "./Data_Availability/renal_paper.wupdtc_pdx_samples.xlsx", sheet = "DNA_RNA")
dnaonly_samples_df <- read_excel(path = "./Data_Availability/renal_paper.wupdtc_pdx_samples.xlsx", sheet = "DNA_only")

# match DNA&RNA samples ---------------------------------------------------
dnarna_samples_count_df <- sampleinfo_df %>%
  filter(SampleID.AcrossDataType %in% dnarna_samples_df$Sample) %>%
  select(SampleID.AcrossDataType) %>%
  table() %>%
  as.data.frame() %>%
  rename(SampleID = '.')
dnarna_samples_df <- merge(x = dnarna_samples_df, y = dnarna_samples_count_df, by.x = "Sample", by.y = "SampleID", all.x = T)

# match DNA samples ---------------------------------------------------
dnaonly_samples_count_df <- sampleinfo_df %>%
  filter(SampleID.AcrossDataType %in% dnaonly_samples_df$Sample) %>%
  select(SampleID.AcrossDataType) %>%
  table() %>%
  as.data.frame() %>%
  rename(SampleID = '.')
dnaonly_samples_df <- merge(x = dnaonly_samples_df, y = dnaonly_samples_count_df, by.x = "Sample", by.y = "SampleID", all.x = T)
