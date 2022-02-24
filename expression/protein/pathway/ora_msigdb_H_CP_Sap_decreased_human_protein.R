# Yige Wu @ WashU 2022 Jan

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input genes
genes_process_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/pathway/map_human_proteins_treated_vs_control_entrezids/20220209.v1/Wilcox.Paired.Each_Treated_Group_vs_CT.20220209.v1.tsv")
## input genes-to-entrez id mapping
genesymbol2entrezid_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/pathway/map_human_proteins_treated_vs_control_entrezids/20220209.v1/genesymbol2entrezid.tsv")
## input MsigDB 
wp2gene1 <- read.gmt("./Resources/Knowledge/Databases/MSigDB/h.all.v7.4.entrez.gmt")
wp2gene2 <- read.gmt("./Resources/Knowledge/Databases/MSigDB/c2.cp.v7.4.entrez.gmt")
wp2gene <- rbind(wp2gene1, wp2gene2)

# get entrez ids ----------------------------------------------------------
genes_background_df <- genes_process_df %>%
  filter(group1 == "Sap") %>%
  mutate(genesymbol = PG.Gene) %>%
  mutate(is_mousegene = str_detect(genesymbol,"[[:lower:]]")) %>%
  filter(!is_mousegene)
## map for markers to be tested
genes_background_df$entrezgene_id <- mapvalues(x = genes_background_df$genesymbol, from = genesymbol2entrezid_df$hgnc_symbol, to = as.vector(genesymbol2entrezid_df$entrezgene_id))
genes_background_df$entrezgene_id <- as.character(genes_background_df$entrezgene_id)
genes_background_df <- genes_background_df %>%
  filter(entrezgene_id != genesymbol) %>%
  filter(!grepl(pattern = "MT\\-", x = genesymbol))
genes_test_df <- genes_background_df %>%
  filter(pvalue < 0.05) %>%
  filter(diff_estimate < 0)
## prepare inputs for clusterprofiler
entrezgene_ids_test <- unique(genes_test_df$entrezgene_id)
entrezgene_ids_universe <- unique(genes_background_df$entrezgene_id)

# test over-representation analysis and gene set enrichment using wikipathway ------------------------------
enricher_out <- tryCatch(expr = enricher(gene = entrezgene_ids_test, TERM2GENE = wp2gene, pvalueCutoff = 1, universe = entrezgene_ids_universe),
                         error = function(e) {warning("ORA failed.");return(NULL)})

if (length(enricher_out) > 0 ) {
  ## convert entrez id to gene symbol
  enricher_out <- setReadable(enricher_out, org.Hs.eg.db, keyType = "ENTREZID")
  enricher_out_all_df <- enricher_out@result
}
# enricher_out_pairwise <- enrichplot::pairwise_termsim(enricher_out)

# plot --------------------------------------------------------------------
# p <- emapplot(x = enricher_out,showCategory = min(50, nrow(enricher_out_all_df[enricher_out_all_df$pvalue < 0.05,])))
# file2write <- paste(dir_out, "emapplot.pdf")
# pdf(file2write, width = 15, height = 10, useDingbats = F)
# print(p)
# dev.off()

# p <- dotplot(object = enricher_out, showCategory=min(50, nrow(enricher_out_all_df[enricher_out_all_df$pvalue < 0.05,])))
# file2write <- paste(dir_out, "dotplot.pdf")
# pdf(file2write, width = 7, height = 4, useDingbats = F)
# print(p)
# dev.off()

# save output -------------------------------------------------------------
# store results
file2write <- paste0(dir_out, "ORA_Results", ".RDS")
saveRDS(object = enricher_out, file = file2write, compress = T)
## store results
file2write <- paste0(dir_out, "ORA_Results", ".tsv")
write.table(x = enricher_out_all_df, file = file2write, quote = F, row.names = F, sep = "\t")

