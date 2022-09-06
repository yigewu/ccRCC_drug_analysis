# Yige Wu @ WashU 2022 Jan

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
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
genes_process_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/rna/test/ttest_paired_diff_rna_1month_combo_vs_single/20220322.v1/mRNA.Ttest.Paired.1month.Combo_vs_single.20220322.v1.tsv")

# get entrez ids ----------------------------------------------------------
genes2convert <- unique(genes_process_df$Name)
## retrieve entrezgene_id
## get biomart
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset(dataset = "hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
genesymbol2entrezid_df <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'),
                                filters = 'hgnc_symbol',
                                values = genes2convert,
                                mart = ensembl)

# write output ------------------------------------------------------------
## store results
file2write <- paste0(dir_out, "genesymbol2entrezid", ".tsv")
write.table(x = genesymbol2entrezid_df, file = file2write, quote = F, row.names = F, sep = "\t")
