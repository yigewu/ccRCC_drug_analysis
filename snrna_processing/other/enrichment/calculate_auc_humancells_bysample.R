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
## input path to the gene rankings
rankings_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/other/enrichment/build_ranking_humancells_bysample/20210216.v1/HumanCells.cell_rankings.paths.20210216.v1.tsv")
## input gene set
geneSets_hallmark <- getGmt("./Resources/Knowledge/Databases/MSigDB/h.all.v7.2.symbols.gmt")
geneSets_kegg <- getGmt("./Resources/Knowledge/Databases/MSigDB/c2.cp.kegg.v7.2.symbols.gmt")
geneSets_react <- getGmt("./Resources/Knowledge/Databases/MSigDB/c2.cp.reactome.v7.2.symbols.gmt")
geneSets_wiki <- getGmt("./Resources/Knowledge/Databases/MSigDB/c2.cp.wikipathways.v7.2.symbols.gmt")

# preprocess --------------------------------------------------------------
# geneSets <- GeneSetCollection(c(geneSets_hallmark, geneSets_kegg, geneSets_react, geneSets_wiki))
geneSets <- geneSets_hallmark

# Calculate enrichment for the gene signatures (AUC) ----------------------
paths_auc <- NULL
for (sampleid_tmp in rankings_df$sampleid) {
  path_rds <- rankings_df$path_rdata[rankings_df$sampleid == sampleid_tmp]
  ## input rds
  cells_rankings <- readRDS(file = path_rds)
  
  ## run AUCell_calcAUC
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, nCores = 4)
  file2write <- paste0(dir_out, sampleid_tmp, ".AUC.rds")
 
  ## store path
  paths_auc <- c(paths_auc, file2write)
}

# write paths -------------------------------------------------------------
paths_df <- data.frame(sampleid = rankings_df$sampleid, path_rankings = rankings_df$path_rdata, path_auc = paths_auc)
table2write <- paste0(dir_out, "HumanCells.AUC.paths.", run_id, ".tsv")
write.table(x = paths_df, file = table2write, quote = F, sep = "\t", row.names = F)

