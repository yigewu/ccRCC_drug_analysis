# Yige Wu @ WashU 2020 Feb
## annotate PDX mutation status within SMGs

# set up libraries and output directory -----------------------------------
## set up working directory and source functions and load libraries
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/ccRCC_Drug_shared.R")
## set run id 
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies --------------------------------------------
## specify the genes to look for mutations
genes <- SMGs[["CCRCC"]]
## input mutation table from the data freeze
maf_datafreeze_tab <- fread("./PDX-Pilot/DataFreeze/SomaticMut/v2.20200126/datafreeze.somaticMut.all.non-silent.meta3.plus.renamedNormal.remMutNearIndel.remPDXmut_extra.remSpeTOmut.tsv", data.table = F)
## input mutation table from batch 8
maf_batch8_tab <- fread("./Ding_Lab/Projects_Current/PDX-WashU/batch8/somaitcMut/pdx/b8.tumorNormal.sm.non-silent.plus.meta2.filteredPDXmut.tsv", data.table = F)

# filter mutation table ----------------------------------------------------
## merge mutations
maf_tab <- rbind(maf_batch8_tab %>%
                   select(Tumor_Sample_Barcode, Hugo_Symbol, HGVSp_Short,  t_alt_count, t_ref_count), 
                 maf_datafreeze_tab %>%
                   select(Tumor_Sample_Barcode, Hugo_Symbol, HGVSp_Short, t_alt_count, t_ref_count))

maf_tab <- maf_tab %>%
  filter(grepl(pattern = "RESL", x = Tumor_Sample_Barcode)) %>%
  unique()
maf_tab$Tumor_Sample_Barcode %>% unique()

# write mutation short amino acid change matrix ---------------------------------------
mut_mat <- get_somatic_mutation_detailed_matrix(pair_tab = genes, maf = maf_tab)
write.table(x = mut_mat, file = paste0(dir_out, "RCC_PDX.Mutation_Matrix.", run_id, ".tsv"), quote = F, row.names = F, sep = "\t")

# write mutation VAF matrix ---------------------------------------
vaf_mat <- get_somatic_mutation_vaf_matrix(pair_tab = genes, maf = maf_tab)
write.table(x = vaf_mat, file = paste0(dir_out, "RCC_PDX.Mutation_VAF_Matrix.", run_id, ".tsv"), quote = F, row.names = F, sep = "\t")

