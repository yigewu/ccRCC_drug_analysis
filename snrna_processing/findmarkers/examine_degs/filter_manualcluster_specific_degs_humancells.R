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
deg_bymodel_bycluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findallmarkers_bymodel_bymanualcluster_humanrcells_8sample_integration_on_katmai/20210216.v1/FindAllMarkers.Wilcox.ByModel.ByCluster.20210216.v1.tsv")

# only keep DEGs that are consistent between RESL5 and RESL10 in each cluster --------
## merge
deg_merged_df <- merge(x = deg_bymodel_bycluster_df %>%
                         filter(model == "RESL5"),
                       y = deg_bymodel_bycluster_df %>%
                         filter(model == "RESL10"), by = c("deg_gene_symbol", "manual_cluster"))
nrow(deg_merged_df)
deg_merged_sig_df <- deg_merged_df %>%
  filter(p_val_adj.x < 0.05 & p_val_adj.y < 0.05) %>%
  filter((avg_logFC.x > 0 & avg_logFC.y > 0) | (avg_logFC.x < 0 & avg_logFC.y < 0)) %>%
  arrange(manual_cluster, p_val_adj.x, desc(avg_logFC.x))
nrow(deg_merged_sig_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ManualCluster_Specific_DEGs.", run_id, ".tsv")
write.table(x = deg_merged_sig_df, file = file2write, quote = F, sep = "\t", row.names = F)
