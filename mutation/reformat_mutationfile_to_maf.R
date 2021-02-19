# Yige Wu @ WashU 2021 Jan

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
mut_df <- fread(data.table = F, input = "./Data_Freeze/v1.dataFreeze.washU_rcc/1.somaticMut/rcc.somaticMut.meta3.20200812.tsv")
## input an example maf file
maf_eg_df <- fread(data.table = F, input = "../ccRCC_snRNA/Resources/Bulk_Processed_Data/Somatic_Variants/ccrcc.somatic.consensus.gdc.umichigan.wu.112918.maf")
colnames(maf_eg_df)

# format ------------------------------------------------------------------
## make sure to have the following columns and the right order
# * temp[0] "Hugo_Symbol"
# * temp[36] "HGVSp_Short"
# * temp[4] "Chromosome"
# * Temp[5] "Start_Position"
# * Temp[10[ "Reference_Allele"
# * Temp[12] "Tumor_Seq_Allele2"
maf_df <- matrix(data = NA, nrow = nrow(mut_df), ncol = ncol(maf_eg_df)) %>% as.data.frame()
colnames(maf_df) <- colnames(maf_eg_df)
colnames_shared <- intersect(colnames(mut_df), colnames(maf_eg_df))
for (colname_tmp in colnames_shared) {
  maf_df[,colname_tmp] <- mut_df[,colname_tmp]
}
maf_df$Start_Position <- mut_df$Start_position
maf_df[1:5, c(1, 5, 6, 11, 13, 37)]

# write output ------------------------------------------------------------
write.table(x = maf_df, file = paste0(dir_out, "RCC_PDX.", run_id, ".maf"), quote = F, row.names = F, sep = "\t")
