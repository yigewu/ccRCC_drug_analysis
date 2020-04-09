# Yige Wu @ WashU 2019 Nov
## for making a table for estiamting VAF from the maf file

# source ------------------------------------------------------------------
setwd(dir = "/Users/yigewu/Box/")
source(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/ccRCC_Drug_shared.R")

# set run id ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input maf files from batch 1-6----------------------------------------------------------
maf_tumor_sample_barcodes <- c("PDX_WUR_014_T", "PDX_WUR_016_T", "PDX_WUR_065_3548_T", "PDX_WUR_065_3549_T")
maf_tumor_sample_ids <- c("RESL5B_4495m", "RESL3C_3317m", "RESL4C_3548m", "RESL4C_3549m")

maf_file1 <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/Resources/PDX_Summary/summary_results/wxs.somaticMut/somaticMut_merged.20190714/u54.wxs_somaticMut_merged.maf.rc.caller", data.table = F)
maf_file1_filtered <- maf_file1 %>%
  filter(Tumor_Sample_Barcode %in% maf_tumor_sample_barcodes) %>%
  mutate(VAF=t_alt_count/t_depth) %>%
  mutate(Case = paste0(str_split_fixed(string = Tumor_Sample_Barcode, pattern = "_", n = 4)[,2], "-", str_split_fixed(string = Tumor_Sample_Barcode, pattern = "_", n = 4)[,3])) %>%
  mutate(Sample = mapvalues(x = Tumor_Sample_Barcode, from = maf_tumor_sample_barcodes, to = maf_tumor_sample_ids)) %>%
  select(Case, Sample, Hugo_Symbol, HGVSp_Short, VAF, Variant_Classification) %>%
  filter(Hugo_Symbol %in% SMGs[["CCRCC"]])

# input maf files from batch 7----------------------------------------------------------
maf_file2 <- fread("./Ding_Lab/Projects_Current/PDX-WashU/batch7/somaticMut/b7.somaticMut.filtered.non-silent.meta2.filteredPDXmut.tsv", data.table = F)
maf_file2_filtered <- maf_file2 %>%
  filter(grepl(x = Tumor_Sample_Barcode, pattern = "RESL")) %>%
  mutate(VAF=t_alt_count/t_depth) %>%
  mutate(Case = paste0(str_split_fixed(string = Tumor_Sample_Barcode, pattern = "-", n = 4)[,2], "-", str_split_fixed(string = Tumor_Sample_Barcode, pattern = "-", n = 4)[,3])) %>%
  mutate(Sample = gsub(x = str_split_fixed(string = Tumor_Sample_Barcode, pattern = "-", n = 4)[,4], pattern = "DRESL", replacement = "RESL")) %>%
  select(Case, Sample, Hugo_Symbol, HGVSp_Short, VAF, Variant_Classification) %>%
  filter(Hugo_Symbol %in% SMGs[["CCRCC"]])

# merge maf files ---------------------------------------------------------
maf_filtered <- rbind(maf_file1_filtered, maf_file2_filtered)

# write table -------------------------------------------------------------
write.table(x = maf_filtered, file = paste0(dir_out, "ccRCC_PDX_SMG_VAF.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)

