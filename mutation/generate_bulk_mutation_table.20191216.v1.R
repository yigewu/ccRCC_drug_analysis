# Yige Wu @ WashU 2019 Dec
## annotate PDX mutation status within SMGs

# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/ccRCC_Drug_shared.R")

# functions ---------------------------------------------------------------
get_somatic_mutation_vaf_matrix <- function(pair_tab, maf) {
  genes4mat <- unique(unlist(pair_tab))
  length(genes4mat)
  
  maf <- maf[maf$Hugo_Symbol %in% genes4mat,]
  nrow(maf)
  maf$sampID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
  maf$vaf <- maf$t_alt_count/(maf$t_alt_count + maf$t_ref_count)
  
  mut_mat <- reshape2::dcast(data = maf, Hugo_Symbol ~ sampID, fun =  function(x) {
    VAF <- paste0(unique(x), collapse = ",")
    return(VAF)
  }, value.var = "vaf", drop=FALSE)
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  return(mut_mat)
}

get_somatic_mutation_detailed_matrix <- function(pair_tab, maf) {
  genes4mat <- unique(unlist(pair_tab))
  length(genes4mat)
  
  maf <- maf[maf$Hugo_Symbol %in% genes4mat,]
  nrow(maf)
  maf$sampID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
  
  mut_mat <- reshape2::dcast(data = maf, Hugo_Symbol ~ sampID, fun =  function(x) {
    variant_class <- paste0(unique(x), collapse = ",")
    return(variant_class)
  }, value.var = "HGVSp_Short", drop=FALSE)
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  return(mut_mat)
}

# set run id --------------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)

# set output directory ----------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# specify the genes to look for mutations ---------------------------------
genes <- SMGs[["CCRCC"]]

# input batch 8 mutation calls --------------------------------------------
mut_batch8_tab <- fread("./Ding_Lab/Projects_Current/PDX-WashU/batch8/somaitcMut/pdx/b8.tumorNormal.sm.non-silent.plus.meta2.filteredPDXmut.tsv", data.table = F)
mut_batch8_tab <- mut_batch8_tab %>%
  filter(grepl(pattern = "RESL", x = Tumor_Sample_Barcode))

# merge mutation table ----------------------------------------------------
maf_tab <- mut_batch8_tab

# generate mutation short amino acid change matrix ---------------------------------------
mut_mat <- get_somatic_mutation_detailed_matrix(pair_tab = genes, maf = maf_tab)
vaf_mat <- get_somatic_mutation_vaf_matrix(pair_tab = genes, maf = maf_tab)
write.table(x = mut_mat, file = paste0(dir_out, "RCC_PDX.Mutation_Matrix.", run_id, ".csv"), quote = F, row.names = F, sep = "\t")
write.table(x = vaf_mat, file = paste0(dir_out, "RCC_PDX.Mutation_VAF_Matrix.", run_id, ".csv"), quote = F, row.names = F, sep = "\t")

