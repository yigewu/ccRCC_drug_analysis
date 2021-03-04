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
## input DEGs by RESL5 and 10
deg_bymodel_bycluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_RESL5vs10_mousecells_bycelltypeshorter_on_katmai/20210302.v2/FindMarkers.Wilcox..RESL5_vs_RESL10.ByCelltypeshorter.tsv")

# make model-specific DEGs in controls --------
deg_ec_df <- deg_bymodel_bycluster_df %>%
  filter(manual_cluster == "Endothelial cells")
## examine cabo targets genes
deg_ec_cabo_df <- deg_ec_df %>%
  mutate(deg_gene_symbol_up = toupper(deg_gene_symbol)) %>%
  filter(deg_gene_symbol_up %in% genes_rtk_cabo) %>%
  filter(p_val_adj < 0.05)
deg_ec_ct_df <- deg_ec_df %>%
  filter(sampleid_group1 == "RESL5E-14541-CT2") %>%
  filter(p_val_adj < 0.05)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "MouseEndothelialCells.", "CT.Model_Specific_DEGs.", run_id, ".tsv")
write.table(x = deg_ec_ct_df, file = file2write, quote = F, sep = "\t", row.names = F)
