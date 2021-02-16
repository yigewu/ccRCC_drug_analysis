# Yige Wu @WashU Apr 2020

# set up libraries and output directory -----------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set time stamp for log file
timestamp <- paste0(run_id, ".", format(Sys.time(), "%H%M%S"))
## getting the path to the current script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
path_this_script <- thisFile()
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
library(AUCell)
library(GSEABase)
library(doMC)
library(doRNG)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)
## set log file
sink(file = paste0(dir_out, "Log.", timestamp, ".txt"))

# input dependencies ------------------------------------------------------
## input gene set
geneSets <- getGmt("./Resources/Knowledge/Databases/MSigDB/h.all.v7.2.symbols.gmt")
## set directory to raw gene count
dir_expmat <- "./Resources/Analysis_Results/snrna_processing/other/generate_filtered_raw_counts_matrix/20210216.v1/"

# preprocess --------------------------------------------------------------
filenames_expmat <- list.files(path = dir_expmat)
filenames_expmat <- filenames_expmat[grepl(pattern = "RESL", x = filenames_expmat)]

# 1. Build gene-expression rankings for each cell -------------------------
paths_rdata <- NULL
sampleids <- NULL
for (filename_expmat in filenames_expmat) {
  ## get sample id
  sampleid_tmp <- str_split_fixed(string = filename_expmat, pattern = "\\.", n = 2)[1,1]
  cat(paste0(sampleid_tmp, "\n"))
  path_expmat <- paste0(dir_expmat, filename_expmat)
  
  ## input data
  exp_df <- fread(data.table = F, input = path_expmat)
  print(exp_df[1:5, 1:4])
  
  ## transform
  exprMatrix <- as.matrix(exp_df[,-1])
  rownames(exprMatrix) <- exp_df$gene_symbol
  print(exprMatrix[1:5, 1:4])
  
  plot2write <- paste0(dir_out, sampleid_tmp, ".png")
  png(plot2write, width = 1000, height = 800, res = 150)
  cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=TRUE)
  dev.off()
  
  ## write output
  rdata2write <- paste0(dir_out, sampleid_tmp, ".cells_rankings.rds")
  saveRDS(cells_rankings, file = rdata2write, compress = T)
  
  ## store paths
  paths_rdata <- c(paths_rdata, rdata2write)
  sampleids <- c(sampleids, sampleid_tmp)
}

# write paths -------------------------------------------------------------
paths_df <- data.frame(sampleid = sampleids, path_rdata = paths_rdata)
table2write <- paste0(dir_out, "HumanCells.cell_rankings.paths.", run_id, ".tsv")
write.table(x = paths_df, file = table2write, quote = F, sep = "\t", row.names = F)
