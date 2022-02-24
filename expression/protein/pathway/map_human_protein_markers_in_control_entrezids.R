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
markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/test/cor_1month_controlproteins_vs_relativetumorvolumevscontrol_bytreatment/20220126.v1/spearman.1month.controlproteins_vs_relativetumorvolume.20220126.v1.tsv")

# get entrez ids ----------------------------------------------------------
markers_df <- markers_df %>%
  filter(abundance_type == "total_protein") %>%
  mutate(genesymbol = str_split_fixed(string = ID, pattern = "_", n = 2)[,1]) %>%
  mutate(is_mousegene = str_detect(genesymbol,"[[:lower:]]")) %>%
  filter(!is_mousegene)
genes2convert <- unique(markers_df$genesymbol)
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
## store results
file2write <- paste0(dir_out, "genesymbol2entrezid", ".tsv")
write.table(x = genesymbol2entrezid_df, file = file2write, quote = F, row.names = F, sep = "\t")
