# Yige Wu @WashU Jul 2020

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

# input dependencies ------------------------------------------------------
## input KEGG
load("./Resources/Knowledge/Databases/2015-08-01_Gene_Set.RData")

# make data frame ---------------------------------------------------------
## get genes from https://www.genome.jp/kegg-bin/show_pathway?hsa04150, which include both mTORC1 and mTORC2
sapagenes_df <- data.frame(genesymbol = KEGG[["hsa04150\tmTOR signaling pathway"]], genesetname = "mTOR signaling pathway", source = "KEGG-hsa04150", category = "Sapasertinib Related")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Sapanisertib_Genes.", run_id, ".tsv")
write.table(x = sapagenes_df, file = file2write, quote = F, row.names = F, sep = "\t")



