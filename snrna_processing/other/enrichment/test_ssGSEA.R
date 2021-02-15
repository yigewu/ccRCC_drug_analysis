devtools::install_github("ncborcherding/escape")


suppressPackageStartupMessages(library(escape))
# suppressPackageStartupMessages(library(dittoSeq))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(Seurat))
pbmc_small <- get("pbmc_small")
pbmc_small <- suppressMessages(UpdateSeuratObject(pbmc_small))
GS <- getGeneSets(library = "H") ## H = hallmark
GS <- getGeneSets(library = "C2") ## H = hallmark

GS_subset <- GS %>%
  filter()