# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
library(clusterProfiler)
library(biomaRt)
library(org.Hs.eg.db)
library(ggplot2)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## get biomart
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
## input wikipathway 
wp2gene <- read.gmt("./Resources/Knowledge/Databases/Wikipathways/wikipathways-20200510-gmt-Homo_sapiens.gmt")
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
## input degs
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/examine_degs/filter_model_specific_bymanualcluster_degs_humancells/20210217.v1/CT.Model_Specific_DEGs.20210217.v1.tsv")

# convert gene symbol to entrez ids ---------------------------------------
## filter
deg_df <- deg_df %>%
  filter(number_up == 0 &  number_down >=3)
genes2convert <- unique(deg_df$deg_gene_symbol)
## retrieve entrezgene_id
genesymbol2entrezid_df <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), 
                                filters = 'hgnc_symbol', 
                                values = genes2convert, 
                                mart = ensembl)
## add entrez ids to the deg table
deg_df$entrezgene_id <- mapvalues(x = deg_df$deg_gene_symbol, from = genesymbol2entrezid_df$hgnc_symbol, to = as.vector(genesymbol2entrezid_df$entrezgene_id))
## filter markers by entrez id
deg_filtered_df <- deg_df %>%
  filter(!is.na(entrezgene_id)) %>%
  filter(entrezgene_id != deg_gene_symbol) %>%
  mutate(avg_logFC = -(mean_avg_logFC))
nrow(deg_df)

## initiate result 
ovaresult_df <- NULL
ovaresult_list <- list()
## get gene list
de_genes <- deg_filtered_df$avg_logFC
de_genes
names(de_genes) <- deg_filtered_df$entrezgene_id
## sort
de_genes <- sort(de_genes, decreasing = T)
de_genes

# test over-representation analysis and gene set enrichment using wikipathway ------------------------------
ora_up_tmp <- tryCatch(expr = enricher(gene = names(de_genes), TERM2GENE = wpid2gene, TERM2NAME = wpid2name, pvalueCutoff = 0.1),
                       error = function(e) {warning("ORA failed.");return(NULL)})
if (length(ora_up_tmp) > 0 ) {
  ## convert entrez id to gene symbol
  ora_up_tmp <- setReadable(ora_up_tmp, org.Hs.eg.db, keyType = "ENTREZID")
  ## store results
  ora_up_tmp_df <- as.data.frame(ora_up_tmp)
  ovaresult_df <- rbind(ovaresult_df, ora_up_tmp_df)
  ## make plots
  p <- cnetplot(x = ora_up_tmp, foldChange=de_genes)
  p <- p + scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0)
  p <- p + labs(color = "Average Expression\nFold Change")
  file2write <- paste0(dir_out, "RESL10_CT_up", ".png")
  png(filename = file2write, width = 1200, height = 1000, res = 150)
  print(p)
  dev.off()
}

# save outputs ------------------------------------------------------------
## save data frames
file2write <- paste0(dir_out, "ORA.", "Result.", run_id, ".tsv")
write.table(x = ovaresult_df, file = file2write, quote = F, sep = "\t", row.names = F)
## save lists
file2write <- paste0(dir_out, "ORA.", "Result.", run_id, ".RDS")
saveRDS(file = file2write, object = ora_up_tmp, compress = T)



