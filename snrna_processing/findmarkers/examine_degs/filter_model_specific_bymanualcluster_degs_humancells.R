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
deg_bymodel_bycluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_RESL5vs10_humancells_bymanualcluster_on_katmai/20210215.v2/FindMarkers.Wilcox..RESL5_vs_RESL10.ByCluster.tsv")

# make model-specific DEGs in controls --------
deg_ct_df <- deg_bymodel_bycluster_df %>%
  filter(sampleid_group1 == "RESL5E-14541-CT2") %>%
  filter(p_val_adj < 0.05)
deg_ct_allcluster_df <- dcast(data = deg_ct_df, formula = deg_gene_symbol ~ manual_cluster, value.var = "avg_logFC")
deg_ct_allcluster_df$number_up <- rowSums(deg_ct_allcluster_df[, 2:5] > 0, na.rm = T)
deg_ct_allcluster_df$number_down <- rowSums(deg_ct_allcluster_df[, 2:5] < 0, na.rm = T)
deg_ct_allcluster_df$mean_avg_logFC <- rowMeans(deg_ct_allcluster_df[, 2:5], na.rm = T)
deg_ct_allcluster_df <- deg_ct_allcluster_df %>%
  arrange(desc(abs(number_up - number_down)), desc(number_down), mean_avg_logFC)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "CT.Model_Specific_DEGs.", run_id, ".tsv")
write.table(x = deg_ct_allcluster_df, file = file2write, quote = F, sep = "\t", row.names = F)
