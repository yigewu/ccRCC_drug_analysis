# Yige Wu @ WashU 2020 Feb

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
## input mutation table from the data freeze
maf_batch1_8_df <- fread("./Resources/Bulk_Processed_Data/Data_Files/batch1_8/SomaticMut/v5.20200408/datafreeze.somaticMut.all.non-silent.meta3.plus.renamedNormal.remMutNearIndel.remPDXmut_extra.remSpeTOmut.tsv.freezedSampleSet.v5.tsv", data.table = F)
# ## input mutation table from the data freeze
# maf_batch7_df <- fread("~/Box/Ding_Lab/Projects_Current/PDX-WashU/batch7/somaticMut/b7.somaticMut.filtered.non-silent.meta2.filteredPDXmut.tsv", data.table = F)
# ## input mutation table from the data freeze
# maf_batch8_df <- fread("~/Box/Ding_Lab/Projects_Current/PDX-WashU/batch8/somaitcMut/pdx/b8.tumorNormal.sm.non-silent.plus.meta2.filteredPDXmut.tsv", data.table = F)
## input mutation table from batch 8
maf_batch9_df <- fread("./Resources/Bulk_Processed_Data/Data_Files/batch9/somaticMut/b9.somaticMut.merged.remMutNearIndel.remPDXmut.meta3.tsv", data.table = F)
## input all ccRCC samples data info
data_status_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_status/write_bulk_data_status/20200526.v1/RCC_PDX_Related_Samples.Batch1_10.Data_Status.20200526.v1.tsv")

# function ----------------------------------------------------------------
get_somatic_mutation_aachange_vaf_matrix <- function(pair_tab, maf) {
  genes4mat <- unique(unlist(pair_tab))
  length(genes4mat)
  
  maf <- maf[maf$Hugo_Symbol %in% genes4mat & maf$Variant_Classification != "Silent",]
  nrow(maf)
  maf$sampID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
  maf$vaf <- maf$t_alt_count/(maf$t_alt_count + maf$t_ref_count)
  maf$HGVSp_sim <- gsub(x = maf$HGVSp_Short, pattern = "p\\.", replacement = "")
  
  maf$aachange_vaf <- paste0(maf$HGVSp_sim, "(", signif(x = maf$vaf, digits = 2), ")")
  mut_mat <- reshape2::dcast(data = maf, Hugo_Symbol ~ sampID, fun =  function(x) {
    aahange_vaf <- paste0(unique(x), collapse = ",")
    return(aahange_vaf)
  }, value.var = "aachange_vaf", drop=FALSE)
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  return(mut_mat)
}

# harmonize the variant classification ------------------------------------
maf_batch1_8_df <- maf_batch1_8_df %>%
  mutate(Variant_Classification.old = Variant_Classification) %>%
  mutate(Variant_Classification = mapvalues(x = Variant_Classification.old, 
                                            from = c("Missense", "Nonsense", "Splice_Site", "Frame_Shift_Del", "In_Frame_Del", "Frame_Shift_Ins", "Nonstop", "In_Frame_Ins"), 
                                            to = c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "Frame_Shift_Del", "In_Frame_Del", "Frame_Shift_Ins", "Nonstop_Mutation", "In_Frame_Ins")))

# filter mutation table ----------------------------------------------------
## get common column names
colnames_mafcommon <- intersect(colnames(maf_batch1_8_df), colnames(maf_batch9_df))
colnames_mafcommon <- intersect(colnames(maf_batch8_df), colnames_mafcommon)
colnames_mafcommon <- colnames_mafcommon[!(colnames_mafcommon %in% c("NCBI_Build", "Mutation_Group"))]
## merge mutations
## filter analysis ids
analysis_ids <- data_status_df$Analysis_ID.WES[!is.na(data_status_df$Analysis_ID.WES)]
maf_merged_df <- rbind(maf_batch9_df[,colnames_mafcommon], 
                       maf_batch1_8_df[,colnames_mafcommon])
maf_merged_df <- maf_merged_df %>%
  filter(Tumor_Sample_Barcode %in% data_status_df$Analysis_ID.WES) %>%
  unique()
## add SampleID.AcrossDataType
maf_merged_df$SampleID.AcrossDataType <- mapvalues(x = maf_merged_df$Tumor_Sample_Barcode, from = data_status_df$Analysis_ID.WES, to = data_status_df$SampleID.AcrossDataType)
## check all samples
maf_merged_df$SampleID.AcrossDataType %>% unique()
## change the tumor barcode
maf_merged_df <- maf_merged_df %>%
  rename(Analysis_ID = Tumor_Sample_Barcode) %>%
  mutate(Tumor_Sample_Barcode = SampleID.AcrossDataType)

# write long data frame ---------------------------------------------------
file2write <- paste0(dir_out, "RCC_PDX.", "MAF.", run_id, ".tsv")
write.table(x = maf_merged_df, file = file2write, quote = F, row.names = F, sep = "\t")

# write mutation short amino acid change matrix ---------------------------------------
mut_mat <- get_somatic_mutation_aachange_vaf_matrix(pair_tab = ccRCC_SMGs, maf = maf_merged_df)
mut_df <- t(mut_mat[,-1]) %>% as.data.frame()
## order
genesymbols_mut <- colnames(mut_df)
mut_df$Case <- rownames(mut_df)
genesymbols_mut_driver <- ccRCC_drivers[ccRCC_drivers %in% genesymbols_mut]
genesymbols_mut_othersmg <- genesymbols_mut[!(genesymbols_mut %in% genesymbols_mut_driver)]
mut_df <- mut_df[, c("Case", genesymbols_mut_driver, genesymbols_mut_othersmg)]
write.table(x = mut_df, file = paste0(dir_out, "RCC_PDX.AAChange_VAF.", run_id, ".tsv"), quote = F, row.names = F, sep = "\t")

# write mutation short amino acid change matrix ---------------------------------------
mut_mat <- get_somatic_mutation_detailed_matrix(pair_tab = ccRCC_SMGs, maf = maf_merged_df)
write.table(x = mut_mat, file = paste0(dir_out, "RCC_PDX.AAChange.", run_id, ".tsv"), quote = F, row.names = F, sep = "\t")

# write mutation VAF matrix ---------------------------------------
vaf_mat <- get_somatic_mutation_vaf_matrix(pair_tab = ccRCC_SMGs, maf = maf_merged_df)
write.table(x = vaf_mat, file = paste0(dir_out, "RCC_PDX.Mutation_VAF_Matrix.", run_id, ".tsv"), quote = F, row.names = F, sep = "\t")

