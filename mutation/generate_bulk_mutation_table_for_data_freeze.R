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
maf_df <- fread(data.table = F, input = "./Data_Freeze/v1.dataFreeze.washU_rcc/1.somaticMut/rcc.somaticMut.meta3.20200812.tsv")

# make mutation short amino acid change matrix ---------------------------------------
mut_mat <- get_somatic_mutation_aachange_vaf_matrix(pair_tab = c(ccRCC_SMGs, "NF2"), maf = maf_df)
mut_df <- t(mut_mat[,-1]) %>% as.data.frame()
## order
genesymbols_mut <- colnames(mut_df)
mut_df$Case <- rownames(mut_df)
genesymbols_mut_driver <- ccRCC_drivers[ccRCC_drivers %in% genesymbols_mut]
genesymbols_mut_othersmg <- genesymbols_mut[!(genesymbols_mut %in% genesymbols_mut_driver)]
mut_df <- mut_df[, c("Case", genesymbols_mut_driver, genesymbols_mut_othersmg)]
write.table(x = mut_df, file = paste0(dir_out, "RCC_PDX.AAChange_VAF.", run_id, ".tsv"), quote = F, row.names = F, sep = "\t")
