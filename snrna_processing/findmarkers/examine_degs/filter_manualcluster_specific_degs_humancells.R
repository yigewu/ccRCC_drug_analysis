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
## input DEGs by manual cluster
deg_bycluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findallmarkers_bymanualcluster_humanrcells_8sample_integration_on_katmai/20210212.v1/FindAllMarkers.Wilcox.ByCluster.20210212.v1.tsv")
## input DEGs by RESL5 and 10
deg_bymodel_bycluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_RESL5vs10_humancells_bymanualcluster_on_katmai/20210215.v1/FindMarkers.Wilcox..RESL5_vs_RESL10.ByCluster.tsv")

# substract the DEGs that are different between RESL5 and RESL10 in each cluster --------
## merge
deg_merged_df <- merge(x = deg_bycluster_df,
                       y = deg_bymodel_bycluster_df)

