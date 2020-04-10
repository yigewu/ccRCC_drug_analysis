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
dir_cellranger_out <- "./Resources/snRNA_Processed_Data/Cell_Ranger/outputs/GRCh38-3.0.0.premrna_and_mm10_premrna/"
## get the sample ids to process
ids_sample <- list.files(path = dir_cellranger_out)
ids_sample

# extract info from metrics_summary.csv -----------------------------------
cellrangermetrics_df <- NULL
for (id_sample in ids_sample) {
  path_cellrangermetrics <- paste0(dir_cellranger_out, id_sample, "/outs/metrics_summary.csv")
  cellrangermetrics_tmp <- fread(input = path_cellrangermetrics, data.table = F)
  cellrangermetrics_tmp <- cellrangermetrics_tmp %>%
    mutate(id_sample = id_sample) %>%
    select(id_sample, `Estimated Number of Cells`,
           `GRCh38-3.0.0.premrna Estimated Number of Cell Partitions`, `mm10 premrna Estimated Number of Cell Partitions`,
           `GRCh38-3.0.0.premrna Median Genes per Cell`, `mm10 premrna Median Genes per Cell`)
  cellrangermetrics_df <- rbind(cellrangermetrics_df, cellrangermetrics_tmp)
}

# write output ------------------------------------------------------------
write.table(x = cellrangermetrics_df, file = paste0(dir_out, "cellrangerfiltered_quality_metrics.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)

