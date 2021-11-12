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
markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/rna/get_foldchange/filter_Cabo_sensitive_resistant_genes_overlapping_proteins/20211110.v2/Cabo_related_genes.20211110.v2.tsv")
background_genes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/rna/get_foldchange/filter_Cabo_sensitive_resistant_genes_overlapping_proteins/20211111.v1/Cabo_resistant_background_genes.20211111.v1.tsv")
## input wikipathway 
wp2gene1 <- read.gmt("./Resources/Knowledge/Databases/MSigDB/h.all.v7.4.entrez.gmt")
wp2gene2 <- read.gmt("./Resources/Knowledge/Databases/MSigDB/c2.cp.v7.4.entrez.gmt")
wp2gene <- rbind(wp2gene1, wp2gene2)

# get entrez ids ----------------------------------------------------------
genes2convert <- unique(background_genes_df$genesymbol)
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

## map for markers to be tested
markers_df$entrezgene_id <- mapvalues(x = markers_df$genesymbol, from = genesymbol2entrezid_df$hgnc_symbol, to = as.vector(genesymbol2entrezid_df$entrezgene_id))
genes_process_df <- markers_df %>%
  filter(gene_category == "Cabo resistant") %>%
  filter(!is.na(entrezgene_id)) %>%
  filter(entrezgene_id != genesymbol) %>%
  filter(!grepl(pattern = "MT\\-", x = genesymbol))
## prepare inputs for clusterprofiler
entrezgene_ids_test <- unique(genes_process_df$entrezgene_id)
entrezgene_ids_universe <- unique(genesymbol2entrezid_df$entrezgene_id); entrezgene_ids_universe <- as.character(entrezgene_ids_universe)

# test over-representation analysis and gene set enrichment using wikipathway ------------------------------
enricher_out <- tryCatch(expr = enricher(gene = entrezgene_ids_test, TERM2GENE = wp2gene, pvalueCutoff = 1, universe = entrezgene_ids_universe),
                         error = function(e) {warning("ORA failed.");return(NULL)})

if (length(enricher_out) > 0 ) {
  ## convert entrez id to gene symbol
  enricher_out <- setReadable(enricher_out, org.Hs.eg.db, keyType = "ENTREZID")
  enricher_out_all_df <- enricher_out@result
}
enricher_out_pairwise <- enrichplot::pairwise_termsim(enricher_out)

# plot --------------------------------------------------------------------
# p <- emapplot(x = enricher_out,showCategory = min(50, nrow(enricher_out_all_df[enricher_out_all_df$pvalue < 0.05,]))) 
# file2write <- paste(dir_out, "emapplot.pdf")
# pdf(file2write, width = 15, height = 10, useDingbats = F)
# print(p)
# dev.off()

p <- dotplot(object = enricher_out, showCategory=min(50, nrow(enricher_out_all_df[enricher_out_all_df$pvalue < 0.01,])))
file2write <- paste(dir_out, "dotplot.pdf")
pdf(file2write, width = 7, height = 4, useDingbats = F)
print(p)
dev.off()


# save output -------------------------------------------------------------
# store results
file2write <- paste0(dir_out, "ORA_Results", ".RDS")
saveRDS(object = enricher_out_pairwise, file = file2write, compress = T)
## store results
file2write <- paste0(dir_out, "ORA_Results", ".tsv")
write.table(x = enricher_out_all_df, file = file2write, quote = F, row.names = F, sep = "\t")

