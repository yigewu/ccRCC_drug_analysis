# Yige Wu @WashU Jul 2020

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

# filter DEGs -------------------------------------------------------------
deg_filtered_df <- deg_df %>%
  filter(pct.1 > 0.1 | pct.2 > 0.1) %>%
  rename(deg_feature_name = deg_gene_symbol) %>%
  mutate(Species_Feature = ifelse(grepl(x = deg_feature_name, pattern = "GRCh38"), "Human", "Mouse")) %>%
  filter(Species_Feature == "Human") %>%
  mutate(deg_gene_symbol = str_split_fixed(string = deg_feature_name, pattern = "GRCh38-3.0.0.premrna-", n = 2)[,2]) %>%
  mutate(Model_ID = str_split_fixed(string = sampleid_group1, pattern = "[A-Z]-", n = 3)[,1]) %>%
  mutate(Treatment_Group1 = str_split_fixed(string = sampleid_group1, pattern = "-", n = 3)[,3]) %>%
  mutate(Treatment_Group1 = gsub(x = Treatment_Group1, pattern = '[1-9]', replacement = "")) %>%
  filter(p_val_adj < 0.05) %>%
  mutate(FoldChange_Direction = ifelse(avg_logFC > 0, "up", "down"))

# overlap DEGs for the same treatment -------------------------------------
deg_filtered_df1 <- deg_filtered_df %>%
  filter(Model_ID == "RESL5") %>%
  select(Model_ID, Treatment_Group1, sampleid_group1, Species_Feature, deg_gene_symbol, FoldChange_Direction, avg_logFC, p_val_adj, pct.1, pct.2, p_val)
deg_filtered_df2 <- deg_filtered_df %>%
  filter(Model_ID == "RESL10") %>%
  select(Model_ID, Treatment_Group1, sampleid_group1, Species_Feature, deg_gene_symbol, FoldChange_Direction, avg_logFC, p_val_adj, pct.1, pct.2, p_val)
## overlap
deg_filtered_overlap_df <- merge(deg_filtered_df1, deg_filtered_df2, by = c("Treatment_Group1", "Species_Feature", "deg_gene_symbol", "FoldChange_Direction"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DEGs_Overalp_Btw_RESL5_vs_RESL10_For_the_Same_Treatment.", run_id, ".tsv")
write.table(x = deg_filtered_overlap_df, file = file2write, sep = "\t", quote = F, row.names = F)

