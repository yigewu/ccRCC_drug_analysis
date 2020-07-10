# Yige Wu @WashU Jun 2020

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
## input DEGs
deg_df1 <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_RESL5_tumorcells_integration_withanchor_on_katmai/20200617.v1/FindMarkers.Wilcox.RESL5.Tumor_cells.integration.withanchor.20200507.v1.Treated_vs_CT..tsv")
deg_df2 <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_RESL10_tumorcells_integration_withanchor_on_katmai/20200702.v1/FindMarkers.Wilcox.RESL10.Tumor_cells.integration.withanchor.20200507.v1.Treated_vs_CT..tsv")
deg_df <- rbind(deg_df1, deg_df2)
## input genes to search for
rccgenes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/write_rcc_genes/20200702.v1/KEGG_RCC.20200702.v1.tsv")
cabogenes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/write_cabozantinib_genes/20200702.v1/Cabozantinib_Genes.20200702.v1.tsv")
sapagenes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/write_sapanisertib_genes/20200702.v1/Sapanisertib_Genes.20200702.v1.tsv")
genes_search_df <- rbind(rccgenes_df, cabogenes_df, sapagenes_df)
genes_search <- unique(genes_search_df$genesymbol)

# filter DEGs -------------------------------------------------------------
deg_filtered_df <- deg_df %>%
  filter(pct.1 > 0.1 | pct.2 > 0.1) %>%
  rename(deg_feature_name = deg_gene_symbol) %>%
  mutate(Species = ifelse(grepl(x = deg_feature_name, pattern = "GRCh38"), "Human", "Mouse")) %>%
  filter(Species == "Human") %>%
  mutate(deg_gene_symbol = str_split_fixed(string = deg_feature_name, pattern = "GRCh38-3.0.0.premrna-", n = 2)[,2]) %>%
  filter(deg_gene_symbol %in% genes_search)


