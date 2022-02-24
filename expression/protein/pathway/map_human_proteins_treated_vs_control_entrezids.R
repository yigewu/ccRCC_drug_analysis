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
genes_process_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/test/wilcox_paired_diff_protein_treated_vs_control/20220209.v1/Wilcox.Paired.Each_Treated_Group_vs_CT.20220209.v1.tsv")

# split by gene -----------------------------------------------------------
idx_rep <- sapply(genes_process_df$PG.Genes, function(string_genes) {
  genes <- str_split(string = string_genes, pattern = ";")[[1]]
  length_genes <- length(genes)
  return(length_genes)
})
genesymbols <- sapply(genes_process_df$PG.Genes, function(string_genes) {
  genes <- str_split(string = string_genes, pattern = ";")[[1]]
  return(genes)
})
genes_process_split_df <- genes_process_df[rep(1:nrow(genes_process_df), idx_rep),]
genes_process_split_df$PG.Gene <- unlist(genesymbols)

# get entrez ids ----------------------------------------------------------
genes_process_filtered_df <- genes_process_split_df %>%
  mutate(is_mousegene = str_detect(PG.Gene,"[[:lower:]]")) %>%
  filter(!is_mousegene)
genes2convert <- unique(genes_process_filtered_df$PG.Gene)
## retrieve entrezgene_id
## get biomart
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
genesymbol2entrezid_df <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'),
                                filters = 'hgnc_symbol',
                                values = genes2convert,
                                mart = ensembl)

# write output ------------------------------------------------------------
## store results
file2write <- paste0(dir_out, "genesymbol2entrezid", ".tsv")
write.table(x = genesymbol2entrezid_df, file = file2write, quote = F, row.names = F, sep = "\t")
file2write <- paste0(dir_out, "Wilcox.Paired.Each_Treated_Group_vs_CT.20220209.v1.tsv")
write.table(x = genes_process_split_df, file = file2write, quote = F, row.names = F, sep = "\t")
