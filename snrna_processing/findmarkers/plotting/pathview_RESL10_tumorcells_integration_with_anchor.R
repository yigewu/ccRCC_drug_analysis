# Yige Wu @WashU Jun 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
source("./ccRCC_drug_analysis/plotting.R")
library(clusterProfiler)
library(pathview)
library(org.Hs.eg.db)
library(dplyr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input DEGs
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_RESL10_tumorcells_integration_withanchor_on_katmai/20200702.v1/FindMarkers.Wilcox.RESL10.Tumor_cells.integration.withanchor.20200507.v1.Treated_vs_CT..tsv")

# map genesymbol to entrez ids --------------------------------------------
## get gene symbol
deg_df <- deg_df %>%
  dplyr::rename(deg_feature_name = deg_gene_symbol) %>%
  dplyr::mutate(Species = ifelse(grepl(x = deg_feature_name, pattern = "GRCh38"), "Human", "Mouse")) %>%
  dplyr::filter(Species == "Human") %>%
  dplyr::mutate(deg_gene_symbol = str_split_fixed(string = deg_feature_name, pattern = "GRCh38-3.0.0.premrna-", n = 2)[,2])
## get unique genes
genes2show <- unique(deg_df$deg_gene_symbol[deg_df$p_val_adj < 0.05])
## map to entrez ids
entrez_ids = bitr(genes2show, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
## map values
deg_df$entrez_id <- mapvalues(x = deg_df$deg_gene_symbol, from = entrez_ids$SYMBOL, to = entrez_ids$ENTREZID)

# plot by sample pair -----------------------------------------------------
for (sampleid_treated in unique(deg_df$sampleid_group1)) {
  dir_out_now <- paste0(dir_out, sampleid_treated, "/")
  dir.create(dir_out_now)
  setwd(dir_out_now)
  ## filter for data to plot
  plot_data_df <- deg_df %>%
    filter(sampleid_group1 == sampleid_treated) %>%
    filter(!is.na(entrez_id))
  ## make named gene list
  geneList <- plot_data_df$avg_logFC
  names(geneList) <- plot_data_df$entrez_id
  for (id_keggpathway in c("hsa05211", "hsa04150", "hsa04370", "hsa04066", "hsa01521", "hsa05200", "hsa04151", "hsa04014")) {
    pathview(gene.data  = geneList,
             pathway.id = id_keggpathway,
             species    = "hsa",
             limit      = list(gene=2, cpd=1))
  }
}

