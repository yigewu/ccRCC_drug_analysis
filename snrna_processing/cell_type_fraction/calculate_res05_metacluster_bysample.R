# Yige Wu @WashU Feb 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
packages = c(
  "rstudioapi",
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
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input data ------------------------------------------------------
barcode2group_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_byres_humancells_8sample_integration_on_katmai/20220905.v1/HumanCells.8sample.Metadata.ByResolution.20220905.v1.tsv")

# calculate fraction ------------------------------------------------------
barcodes_by_cluster_sample_df <- barcode2group_df %>%
  select(integrated_snn_res.0.5, orig.ident) %>%
  rename(Id_Sample = orig.ident) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 0)
barcodes_by_sample_df <- barcode2group_df %>%
  rename(Id_Sample = orig.ident) %>%
  select(Id_Sample) %>%
  table() %>%
  as.data.frame() %>%
  rename(Id_Sample = '.')
## merge number of barcodes
barcodes_by_cluster_sample_df <- merge(barcodes_by_cluster_sample_df, barcodes_by_sample_df, by = c("Id_Sample"), suffixes = c("_Cluster_Sample", "_Sample"))
## calculate fraction
barcodes_by_cluster_sample_df <- barcodes_by_cluster_sample_df %>%
  mutate(Fraction_Cluster_Sample = Freq_Cluster_Sample/Freq_Sample)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "MetaClusterFraction.", run_id, ".tsv")
write.table(x = barcodes_by_cluster_sample_df, file = file2write, sep = "\t", quote = F, row.names = F)

