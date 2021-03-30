# Yige Wu @WashU Mar 2021

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
## input DEGS between RESL10 vs RESL5
deg_bymodel_snRNA_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_RESL5vs10_mousecells_bycelltypeshorter_on_katmai/20210302.v2/FindMarkers.Wilcox..RESL5_vs_RESL10.ByCelltypeshorter.tsv")
## input bulk protein data
deg_bymodel_bulkpro_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/get_RESL10_vs_RESL5_CT_protein_diff/20210304.v1/RESL10_vs_RESL5.Protein.20210304.v1.tsv")

# filter ------------------------------------------------------------------
deg_filtered_df <- deg_bymodel_snRNA_df %>%
  filter(sampleid_group1 == "RESL5E-14541-CT2") %>%
  filter(manual_cluster == "Endothelial cells") %>%
  filter(p_val_adj < 0.05) %>%
  filter(pct.1 > pct.2) %>%
  filter(avg_logFC > 0)
nrow(deg_filtered_df)

# merged ------------------------------------------------------------------
deg_merged_df <- merge(x = deg_filtered_df,
                       y = deg_bymodel_bulkpro_df %>%
                         select(PG.Gene, PG.ProteinName, PG.ProteinAccession, PG.ProteinNames, 
                                RESL10_vs_RESL5.protein.2m_CT, log2Intensity_FC_RESL10vs5_2m_CT, 
                                RESL10_vs_RESL5.protein.1m_CT, log2Intensity_FC_RESL10vs5_1m_CT), 
                       by.x = c("deg_gene_symbol"), 
                       by.y = c("PG.Gene"), all.x = T)
deg_merged_df <- deg_merged_df %>%
  arrange(p_val_adj)
nrow(deg_merged_df[!is.na(deg_merged_df$RESL10_vs_RESL5.protein.2m_CT),])
# filter by bulk protein level --------------------------------------------
deg_pro_filtered_df <- deg_merged_df %>%
  filter(RESL10_vs_RESL5.protein.2m_CT == "Down")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Mouse.Endothelial.RESL5_vs_RESL10_CT.High.", run_id, ".tsv")
write.table(x = deg_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "Mouse.Endothelial.RESL5_vs_RESL10_CT.High.ProteinFiltered.", run_id, ".tsv")
write.table(x = deg_pro_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)


