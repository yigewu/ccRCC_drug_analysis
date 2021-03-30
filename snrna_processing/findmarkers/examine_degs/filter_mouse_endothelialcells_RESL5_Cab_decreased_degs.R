# Yige Wu @WashU Mar 2021
## these DEGs will be down-regulated in RESL5-Cab-treated vs RESL5-CT 

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
## input DEGs between treated vs control
deg_bytreatment_snRNA_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_TreatedvsCT_mousecells_bycelltypeshorter_on_katmai/20210323.v1/FindMarkers.Wilcox.Treated_vs_CT.ByCellType.tsv")
## input bulk protein data
deg_bytreatment_bulkpro_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/get_treated_vs_CT_protein_diff/20210316.v1/Treated_vs_CT.Protein.20210316.v1.tsv")

# get genes that are down-regulated in RESL10-Sap-treated vs RESL10-CT, higher in RESL10-CT than RESL5-CT ------------------------------------------------------------------
deg_cab_dec_resl5_df <- deg_bytreatment_snRNA_df %>%
  filter(sampleid_group1 == "RESL5E-14539-Cabo2") %>%
  filter(manual_cluster == "Endothelial cells") %>%
  filter(avg_logFC < 0) %>%
  filter(p_val_adj < 0.05) %>%
  # filter(pct.1 < pct.2) %>%
  arrange(avg_logFC)
nrow(deg_cab_dec_resl5_df)

# add bulk RNA fold change and bulk protein fold change -------------------
deg_merged_df <- deg_cab_dec_resl5_df
deg_merged_df <- merge(x = deg_merged_df,
                       y = deg_bytreatment_bulkpro_df %>%
                         select(PG.Gene, PG.ProteinName, PG.ProteinAccession, PG.ProteinNames, log2Intensity_FC_RESL5_2m_CabvsCT, log2Intensity_FC_RESL5_1m_CabvsCT), 
                       by.x = c("deg_gene_symbol"), 
                       by.y = c("PG.Gene"), all.x = T)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Mouse.Endothelial.RESL5_Cabo_vs_CT.Decreased.", run_id, ".tsv")
write.table(x = deg_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)

