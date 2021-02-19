# Yige Wu @WashU Feb 2020

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
deg_bymodel_ct_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/examine_degs/filter_model_specific_bymanualcluster_degs_humancells/20210217.v1/CT.Model_Specific_DEGs.20210217.v1.tsv")
## input bulk RNA DEGs
deg_bulkrna_df <- fread(data.table = F, input = "./Analysis/bulk.DEGs/resl10_vs_resl5/resl10_vs_resl5.deg.20210219.out")

# filter RESL10 up DEGs------------------------------------------------------------------
deg_resl10_up_ct_df <- deg_bymodel_ct_df %>%
  filter(number_down >= 3 & number_up == 0)
deg_resl10_up_ct_df <- merge(x = deg_resl10_up_ct_df, 
                             y = deg_bulkrna_df %>%
                               filter(RESL10_vs_RESL5 == "Up") %>%
                               filter(Status %in% c("Control")) %>%
                               filter(Group == "2-month"),
                             by.x = c("deg_gene_symbol"),
                             by.y = c("Gene"),
                             all.x = T)
deg_resl10_up_ct_df <- deg_resl10_up_ct_df %>%
  arrange(mean_avg_logFC)