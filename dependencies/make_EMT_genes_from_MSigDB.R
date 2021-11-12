# Yige Wu @ WashU 2021 Nov

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
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)

# input dependencies ------------------------------------------------------
entrezid2geneset_df <- read.gmt("./Resources/Knowledge/Databases/MSigDB/h.all.v7.4.entrez.gmt")

# get biomart -------------------------------------------------------------
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)

# get entrez ids ----------------------------------------------------------
genes2convert <- unique(entrezid2geneset_df$gene[entrezid2geneset_df$term == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"])
## retrieve entrezgene_id
genesymbol2entrezid_df <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), 
                                filters = 'entrezgene_id', 
                                values = genes2convert, 
                                mart = ensembl)

# save output -------------------------------------------------------------
## store results
file2write <- paste0(dir_out, "MSigDB_EMT_genesymbols", ".tsv")
write.table(x = genesymbol2entrezid_df, file = file2write, quote = F, row.names = F, sep = "\t")

