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
maf_df <- fread(data.table = F, input = "./Data_Freeze/v1.dataFreeze.washU_rcc/1.somaticMut/rcc.somaticMut.meta3.20201007.tsv")

# make mutation short amino acid change matrix ---------------------------------------
genes4mat <- c(ccRCC_SMGs, "NF2")
maf_filtered_df <- maf_df[maf_df$Hugo_Symbol %in% genes4mat & maf_df$Variant_Classification != "Silent",]
maf_filtered_df$vaf <- maf_filtered_df$t_alt_count/(maf_filtered_df$t_alt_count + maf_filtered_df$t_ref_count)
maf_filtered_df$HGVSp_sim <- gsub(x = maf_filtered_df$HGVSp_Short, pattern = "p\\.", replacement = "")
maf_filtered_df$aachange_vaf <- paste0(maf_filtered_df$HGVSp_sim, "(", signif(x = maf_filtered_df$vaf, digits = 2), ")")
mut_mat <- reshape2::dcast(data = maf_filtered_df, Hugo_Symbol ~ Tumor_Sample_Barcode, fun =  function(x) {
  aahange_vaf <- paste0(unique(x), collapse = ",")
  return(aahange_vaf)
}, value.var = "aachange_vaf", drop=FALSE)
rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
mut_df <- t(mut_mat[,-1]) %>% as.data.frame()
## order
genesymbols_mut <- colnames(mut_df)
mut_df$Case <- rownames(mut_df)
genesymbols_mut_driver <- ccRCC_drivers[ccRCC_drivers %in% genesymbols_mut]
genesymbols_mut_othersmg <- genesymbols_mut[!(genesymbols_mut %in% genesymbols_mut_driver)]
mut_df <- mut_df[, c("Case", genesymbols_mut_driver, genesymbols_mut_othersmg)]
write.table(x = mut_df, file = paste0(dir_out, "RCC_PDX.AAChange_VAF.", run_id, ".tsv"), quote = F, row.names = F, sep = "\t")
