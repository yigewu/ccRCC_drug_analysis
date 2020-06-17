# Yige Wu @ WashU 2020 Feb
## annotate PDX CNV status

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies --------------------------------------------
## input all ccRCC samples data info
data_status_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_status/write_bulk_data_status/20200526.v1/RCC_PDX_Related_Samples.Batch1_10.Data_Status.20200526.v1.tsv")
## input the gistic arm-level data
genecnv_df1 <- fread(input = "~/Box/PDX-Pilot/DataFreeze/CopyNumber/v2.cnvkit.seg5kb.v20200114/v2.tumorOnly.cnvkit_gene-level.5kb.all.merged.filtered.tsv", data.table = F)
genecnv_df2 <- fread(input = "~/Box/PDX-Pilot/DataFreeze/CopyNumber/v2.cnvkit.seg5kb.v20200114/v2.tumorNormal.cnvkit_gene-level.5kb.all.merged.filtered.tsv", data.table = F)
genecnv_df3 <- fread(input = "./Resources/Bulk_Processed_Data/Data_Files/batch9/cnv/b9.cnvkit_segment.all.merged.gene-level.tsv", data.table = F)

# filter results ----------------------------------------------------------
## filter by RCC samples
genecnv_df1_filtered <- genecnv_df1 %>%
  filter(sample %in% data_status_df$Analysis_ID.WES)
genecnv_df2_filtered <- genecnv_df2 %>%
  filter(sample %in% data_status_df$Analysis_ID.WES)
genecnv_df3_filtered <- genecnv_df3 %>%
  filter(sample %in% data_status_df$Analysis_ID.WES) %>%
  rename(cn = segment_cn) %>%
  rename(chromosome = chr) %>%
  rename(start = segment_start) %>%
  rename(end = segment_end) %>%
  rename(log2 = segment_log2) %>%
  rename(depth = segment_depth) %>%
  rename(weight = segment_weight)
## get common columns
colnames_common <- intersect(colnames(genecnv_df1_filtered), colnames(genecnv_df3_filtered))
## combine
genecnv_merged_df <- rbind(genecnv_df1_filtered[,colnames_common], 
                           genecnv_df2_filtered[,colnames_common], 
                           genecnv_df3_filtered[,colnames_common])

# write output ------------------------------------------------------------
write.table(x = genecnv_merged_df, file = paste0(dir_out, "RCC_PDX.Gene_Level_CNV.", run_id, ".tsv"), quote = F, row.names = F, sep = "\t")


