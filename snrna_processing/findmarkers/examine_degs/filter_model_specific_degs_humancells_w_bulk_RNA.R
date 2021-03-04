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
deg_snRNA_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/examine_degs/filter_model_specific_bymanualcluster_degs_humancells/20210217.v1/CT.Model_Specific_DEGs.20210217.v1.tsv")
## input bulk RNA DEGs
deg_bulkrna_df <- fread(data.table = F, input = "./Analysis/bulk.DEGs/resl10_vs_resl5/resl10_vs_resl5.model_level.deg.20210219.out")

# filter and merge --------------------------------------------------------
deg_snRNA_df <- deg_snRNA_df %>%
  mutate(RESL10_vs_RESL5 = ifelse(number_down >= 3 & number_up == 0, "Up",
                                  ifelse(number_up >= 3 & number_down == 0, "Down", "Other")))
deg_merged_df <- merge(x = deg_snRNA_df %>%
                         filter(RESL10_vs_RESL5 != "Other"),
                       y = deg_bulkrna_df, by.x = c("deg_gene_symbol", "RESL10_vs_RESL5"), by.y = c("Gene", "RESL10_vs_RESL5"), all.x = T)
deg_merged_df <- deg_merged_df %>%
  arrange(RESL10_vs_RESL5, -log2FC)

nrow(deg_merged_df[!is.na(deg_merged_df$log2FC) & deg_merged_df$RESL10_vs_RESL5 == "Up",])

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "CT.Model_Specific_DEGs.", run_id, ".tsv")
write.table(x = deg_ct_allcluster_df, file = file2write, quote = F, sep = "\t", row.names = F)
