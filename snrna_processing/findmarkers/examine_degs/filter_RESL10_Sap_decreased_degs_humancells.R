# Yige Wu @WashU Mar 2021
## these DEGs will be down-regulated in RESL10-Sap-treated vs RESL10-CT 
## they are expressedin a lower %  in RESL5 control than RESL10 control

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
deg_bytreatment_snRNA_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_TreatedvsCT_humancells_on_katmai/20210304.v1/FindMarkers.Wilcox.Treated_vs_CT.tsv")
## input DEGS between RESL10 vs RESL5
deg_bymodel_snRNA_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_RESL10vs5_humancells_on_katmai/20210304.v2/FindMarkers.Wilcox..RESL10_vs_RESL5.tsv")
## input bulk RNA DEGs
deg_bulkrna_df <- fread(data.table = F, input = "./Analysis/bulk.DEGs/resl10_vs_resl5/resl10_vs_resl5.catalog_level.deg.20210219.out")
## input bulk protein data
deg_bymodel_bulkpro_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/get_RESL10_vs_RESL5_CT_protein_diff/20210304.v1/RESL10_vs_RESL5.Protein.20210304.v1.tsv")
deg_bytreatment_bulkpro_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/get_treated_vs_CT_protein_diff/20210309.v1/Treated_vs_CT.Protein.20210309.v1.tsv")

# get genes that are down-regulated in RESL10-Sap-treated vs RESL10-CT, higher in RESL10-CT than RESL5-CT ------------------------------------------------------------------
deg_sap_dec_resl10_df <- deg_bytreatment_snRNA_df %>%
  filter(sampleid_group1 == "RESL10F-12465-Sap2") %>%
  filter(avg_logFC < 0) %>%
  filter(p_val_adj < 0.05) %>%
  filter(pct.1 < pct.2) %>%
  arrange(avg_logFC)
nrow(deg_sap_dec_resl10_df)
deg_resl10_high_df <- deg_bymodel_snRNA_df %>%
  filter(sampleid_group1 == "RESL10F-12462-CT2" & avg_logFC > 0 & p_val_adj < 0.05)
## filter to those that are higher in RESL10 CT vs RESL5 CT
deg_sap_dec_resl10_filtered_df <- deg_sap_dec_resl10_df %>%
  filter(deg_gene_symbol %in% deg_resl10_high_df$deg_gene_symbol)
nrow(deg_sap_dec_resl10_filtered_df)

# add bulk RNA fold change and bulk protein fold change -------------------
deg_bulkrna_2mct_df <- deg_bulkrna_df %>%
  filter(Group == "2-month") %>%
  filter(Status == "Control")
deg_merged_df <- merge(x = deg_sap_dec_resl10_filtered_df,
                       y = deg_bulkrna_2mct_df %>%
                         mutate(log2FC.bulkRNA = log2FC) %>%
                         mutate(RESL10_vs_RESL5.bulkRNA = RESL10_vs_RESL5) %>%
                         select(Gene, log2FC.bulkRNA, RESL10_vs_RESL5.bulkRNA), 
                       by.x = c("deg_gene_symbol"), 
                       by.y = c("Gene"), all.x = T)
deg_merged_df <- merge(x = deg_merged_df,
                       y = deg_bymodel_bulkpro_df %>%
                         select(PG.Gene, PG.ProteinName, PG.ProteinAccession, PG.ProteinNames, 
                                RESL10_vs_RESL5.protein.2m_CT, log2Intensity_FC_RESL10vs5_2m_CT, 
                                RESL10_vs_RESL5.protein.1m_CT, log2Intensity_FC_RESL10vs5_1m_CT), 
                       by.x = c("deg_gene_symbol"), 
                       by.y = c("PG.Gene"), all.x = T)
deg_merged_df <- merge(x = deg_merged_df,
                       y = deg_bytreatment_bulkpro_df, 
                       by.x = c("deg_gene_symbol", "PG.ProteinName", "PG.ProteinAccession", "PG.ProteinNames"), 
                       by.y = c("PG.Gene", "PG.ProteinName", "PG.ProteinAccession", "PG.ProteinNames"), all.x = T)

deg_merged_df <- deg_merged_df %>%
  arrange(avg_logFC) %>%
  mutate(RESL10_vs_RESL5.protein.2m_CT = ifelse(is.na(RESL10_vs_RESL5.protein.2m_CT), "NA", RESL10_vs_RESL5.protein.2m_CT)) %>%
  mutate(Sap_vs_CT.protein.2m_RESL10 = ifelse(is.na(log2Intensity_FC_RESL10_2m_SapvsCT), "NA", 
                                              ifelse(log2Intensity_FC_RESL10_2m_SapvsCT > 0, "Up", "Down"))) %>%
  mutate(Sap_vs_CT.protein.1m_RESL10 = ifelse(is.na(log2Intensity_FC_RESL10_1m_SapvsCT), "NA", 
                                              ifelse(log2Intensity_FC_RESL10_1m_SapvsCT > 0, "Up", "Down"))) %>%
  mutate(RESL10_vs_RESL5.bulkRNA = ifelse(is.na(RESL10_vs_RESL5.bulkRNA), "NA", RESL10_vs_RESL5.bulkRNA))

deg_merged_filtered_df1 <- deg_merged_df %>%
  filter(RESL10_vs_RESL5.protein.2m_CT != "Down") %>%
  filter(Sap_vs_CT.protein.2m_RESL10 != "Up") %>%
  # filter(RESL10_vs_RESL5.protein.1m_CT != "Down") %>%
  # filter(Sap_vs_CT.protein.1m_RESL10 != "Up") %>%
  filter(pct.1 < pct.2) %>%
  filter(RESL10_vs_RESL5.bulkRNA != "Down")
nrow(deg_merged_filtered_df1)
table(deg_merged_filtered_df1 %>%
        select(Sap_vs_CT.protein.2m_RESL10, RESL10_vs_RESL5.protein.2m_CT, RESL10_vs_RESL5.bulkRNA))

deg_merged_filtered_df2 <- deg_merged_filtered_df1 %>%
  filter(RESL10_vs_RESL5.protein.2m_CT == "Up") %>%
  filter(Sap_vs_CT.protein.2m_RESL10 == "Down") %>%
  arrange(avg_logFC)

deg_merged_filtered_df3 <- deg_merged_filtered_df1 %>%
  filter(RESL10_vs_RESL5.protein.2m_CT == "Up") %>%
  filter(Sap_vs_CT.protein.2m_RESL10 == "Down") %>%
  filter(RESL10_vs_RESL5.bulkRNA == "Up")


# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DEGs_Decreased_in_RESL10_Sap_vs_CT.", run_id, ".tsv")
write.table(x = deg_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)
