# Yige Wu @ WashU 2022 Mar
## refer to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6075733/

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
packages = c(
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
source("./ccRCC_drug_analysis/functions.R")

# input dependencies ------------------------------------------------------
## input sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/0.sample_info/v3.20210116/RCC_PDX_Samples.20210115.v2.xlsx")
## input CNV
cnv_arm_df <- fread(data.table = F, input = "./Data_Freeze/v1.dataFreeze.washU_rcc/2.copyNumber/gistic2.broad_values_by_arm.20200818.txt")

# process -----------------------------------------------------------------
samples_plot_df <- sampleinfo_df %>%
  filter(ModelID %in% paste0("RESL", c(3, 4, 5, 10, 11, 12))) %>%
  # filter(DataType == "WES") %>%
  filter(Group == "PDX") %>%
  filter(ShortTag %in% c("Baseline", "Control", "Treated.Cabo", "Treated.Sap", "Treated.Cabo+Sap")) %>%
  mutate(column_id = ModelID)

## process CNV arm-level data
cnv_byarm_bysample_df <- data.frame(t(cnv_arm_df[,-1])); colnames(cnv_byarm_bysample_df) <- cnv_arm_df$`Chromosome Arm`; cnv_byarm_bysample_df$Analysis_ID = colnames(cnv_arm_df)[-1]
cnv_byarm_bysample_df <- cnv_byarm_bysample_df %>%
  filter(Analysis_ID %in% samples_plot_df$Analysis_ID)
cnv_byarm_bysample_df$ModelID <- mapvalues(x = cnv_byarm_bysample_df$Analysis_ID, from = sampleinfo_df$Analysis_ID, to = as.vector(sampleinfo_df$ModelID))
cnv_byarm_bymodel_df <- cnv_byarm_bysample_df %>%
  group_by(ModelID) %>%
  rename(log2cr.3p = '3p') %>%
  rename(log2cr.5q = '5q') %>%
  rename(log2cr.14q = '14q') %>%
  rename(log2cr.7p = '7p') %>%
  rename(log2cr.7q = '7q') %>%
  rename(log2cr.17p = '17p') %>%
  rename(log2cr.17q = '17q') %>%
  summarise(avg_log2cr.3p = round(mean(log2cr.3p), digits = 2),
            avg_log2cr.5q = round(mean(log2cr.5q), digits = 2),
            avg_log2cr.14q = round(mean(log2cr.14q), digits = 2),
            avg_log2cr.7 = round(mean(log2cr.7p+log2cr.7q), digits = 2),
            avg_log2cr.17 = round(mean(log2cr.17p+log2cr.17q), digits = 2)) %>%
  arrange(factor(ModelID, levels = paste0("RESL", c(3, 4, 5, 10, 11, 12))))

# write output ------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
## write file
file2write <- paste0(dir_out, "Arm_Level.log2CopyRatio.ByModel.", run_id, ".tsv")
write.table(x = cnv_byarm_bymodel_df, file = file2write, quote = F, sep = "\t", row.names = F)

