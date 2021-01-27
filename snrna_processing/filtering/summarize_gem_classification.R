# Yige Wu @WashU Apr 2020
## for making a table for selected quality metrics after cell ranger processing, before filtering

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
library(Seurat)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## set directory for cell ranger output
dir_cellranger_out <- "./Resources/snRNA_Processed_Data/Cell_Ranger/outputs/cellranger-5.0.1_Ref-2020-A/"
## get the sample ids to process
ids_sample <- list.files(path = dir_cellranger_out)
ids_sample <- ids_sample[grepl(pattern = 'RESL', x = ids_sample)]
ids_sample

# extract info from metrics_summary.csv -----------------------------------
gemclass_summary_df <- NULL

for (id_sample_tmp in ids_sample) {
  path_file <- paste0(dir_cellranger_out, id_sample_tmp, "/", id_sample_tmp, "/outs/analysis/gem_classification.csv")
  gemclass_tmp_df <- fread(input = path_file, data.table = F)
  gemclass_summary_tmp_df <- gemclass_tmp_df %>%
    select(call) %>%
    table() %>%
    as.data.frame() %>%
    rename(gem_classification = '.')
  
  gemclass_summary_tmp_2merge_df <- as.data.frame(t(gemclass_summary_tmp_df$Freq))
  colnames(gemclass_summary_tmp_2merge_df) <- gemclass_summary_tmp_df$gem_classification
  gemclass_summary_tmp_2merge_df <- gemclass_summary_tmp_2merge_df %>%
    mutate(id_sample = id_sample_tmp) %>%
    select(id_sample, GRCh38, mm10, Multiplet)
  gemclass_summary_df <- rbind(gemclass_summary_df, gemclass_summary_tmp_2merge_df)
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "cellrangerfiltered_gem_classification_summary.", run_id, ".tsv")
write.table(x = gemclass_summary_df, file = file2write, quote = F, sep = "\t", row.names = F)

