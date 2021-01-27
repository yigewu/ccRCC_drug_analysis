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
cellrangermetrics_df <- NULL
cellrangermetrics_sub_df <- NULL

for (id_sample_tmp in ids_sample) {
  path_cellrangermetrics <- paste0(dir_cellranger_out, id_sample_tmp, "/", id_sample_tmp, "/outs/metrics_summary.csv")
  cellrangermetrics_tmp <- fread(input = path_cellrangermetrics, data.table = F)
  cellrangermetrics_tmp <- cellrangermetrics_tmp %>%
    mutate(id_sample = id_sample_tmp)
    
  cellrangermetrics_df <- rbind(cellrangermetrics_df, cellrangermetrics_tmp)
  
  cellrangermetrics_sub_tmp <- cellrangermetrics_tmp %>%
    select(id_sample, `Estimated Number of Cells`,
           `GRCh38 Estimated Number of Cell Partitions`, `mm10 Estimated Number of Cell Partitions`,
           `GRCh38 Median Genes per Cell`, `mm10 Median Genes per Cell`)
  cellrangermetrics_sub_df <- rbind(cellrangermetrics_sub_df, cellrangermetrics_sub_tmp)
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "cellrangerfiltered_metrics_summary.", run_id, ".tsv")
write.table(x = cellrangermetrics_df, file = file2write, quote = F, sep = "\t", row.names = F)

file2write <- paste0(dir_out, "cellrangerfiltered_metrics_summary.subset.", run_id, ".tsv")
write.table(x = cellrangermetrics_sub_df, file = file2write, quote = F, sep = "\t", row.names = F)

