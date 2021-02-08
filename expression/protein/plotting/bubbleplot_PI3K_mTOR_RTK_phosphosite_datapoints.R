# Yige Wu @ WashU 2021 Jan

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
source("./ccRCC_drug_analysis/plotting.R")
library(ComplexHeatmap)
## set run id
version_tmp <- 3
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the protein data
exp_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/preprocess/filter_and_transform_DIA_phosphorylation_data/20210205.v1/RCC_PDX.DIA_Phosphopeptide.Log2.20210205.v1.tsv")
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")

# set parameters ----------------------------------------------------------
genes_filter <- c(genes_pi3k_mtor, genes_rtk_cabo)

# make plot data ----------------------------------------------------------
## filter by gene
idx_keep <- sapply(exp_df$PG.Genes, function(gene_string, genes_filter_vec) {
  genes_vec <- str_split(string = gene_string, pattern = ";")[[1]]
  idx_keep_human <- any(genes_vec %in% genes_filter_vec)
  idx_keep_mouse <- any(genes_vec %in% tolower(genes_filter_vec))
  idx_keep_tmp <- (idx_keep_human || idx_keep_mouse)
  return(idx_keep_tmp)
}, genes_filter_vec = genes_filter)
plot_data_df <- exp_df[idx_keep,]

# plot --------------------------------------------------------------------


# write outpu -------------------------------------------------------------


