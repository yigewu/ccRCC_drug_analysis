# Yige Wu @WashU Mar 2021
## these DEGs will be down-regulated in RESL10-Sap-treated vs RESL10-CT (preferably similarly in multiple clusters)
## they are expressed in RESL5 control in a lower % than RESL10 control

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
deg_treated_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_TreatedvsCT_humancells_bymanualcluster_on_katmai/20210216.v1/FindMarkers.Wilcox.Treated_vs_CT.ByCluster.tsv")
## input DEGS between RESL10 vs RESL5
deg_bymodel_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/clusterprofiler/run_clusterprofiler_on_degs_humancells_RESL10_vs_RESL5_bulk_rna_pro_filtered/20210303.v1/DEG.Filtered.20210303.v1.tsv")

# filter ------------------------------------------------------------------
## merge by treatement group
deg_bycluster_df <- deg_treated_df %>%
  mutate(treatment_group = str_split_fixed(string = sampleid_group1, pattern = "\\-", n = 3)[,3])
deg_merged_df <- merge(x = deg_treated_df %>%
                         filter(grepl(x = sampleid_group1, pattern = "RESL10")),
                       y = deg_treated_df %>%
                         filter(grepl(x = sampleid_group1, pattern = "RESL5")),
                       by = c("deg_gene_symbol", "manual_cluster", "treatment_group"), suffixes = c(".RESL10", ".RESL5"))
## filter by cluster
deg_filtered_df <- deg_merged_df %>%
  filter(sampleid_group1.RESL10 == "RESL10F-12465-Sap2") %>%
  # filter(manual_cluster == 1) %>%
  filter(manual_cluster %in% c(1, 2, 3, 4, 5)) %>%
  filter(p_val_adj.RESL10 < 0.05 & avg_logFC.RESL10 < 0) %>%
  filter(pct.2.RESL10 > pct.2.RESL5) %>%
  arrange((pct.2.RESL5 - pct.2.RESL10))


intersect(deg_filtered_df$deg_gene_symbol, deg_bymodel_df$deg_gene_symbol[deg_bymodel_df$RESL10_vs_RESL5 == "Up"])
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ManualCluster_Specific_DEGs.", run_id, ".tsv")
write.table(x = deg_merged_sig_df, file = file2write, quote = F, sep = "\t", row.names = F)
