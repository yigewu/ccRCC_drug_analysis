# Yige Wu @WashU Sep 2022
## for making dimplot for RESL10_4sample_integration

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
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_Drug/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "Seurat"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
## set output directory
source("./ccRCC_drug_analysis/functions.R")
dir_out_parent <- makeOutDir_katmai(path_this_script)
dir_out <- paste0(dir_out_parent, run_id, "/")
dir.create(dir_out)

# input data --------------------------------------------------------
## set the path to the rds file for integrated object
path_rds <- "./Resources/Analysis_Results/snrna_processing/integration/run_humancells_8sample_integration_withanchor_on_katmai/20210208.v1/Humancells_8sample_integration.withanchor.20210208.v1.RDS"
## input RDS file
srat <- readRDS(file = path_rds)
print("Finish reading the RDS file!\n")


# plot --------------------------------------------------------------------
srat@meta.data[,"MC_name"] = paste0("MC", as.numeric(srat@meta.data[,"integrated_snn_res.0.5"]))
Idents(srat) = "MC_name"
for (gene in c("ABCA1", "C3", "NAV2", "NEAT1", "PLD1", "ROR1")) {
  pdf(paste0(dir_out, gene, ".pdf"), width = 8, height = 4, useDingbats = F)
  Seurat::VlnPlot(object = srat, features = gene)
  dev.off()
}
