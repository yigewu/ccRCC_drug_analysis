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
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "Seurat",
  "future",
  "future.apply"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
# set up future for parallelization
plan("multiprocess", workers = 6)
options(future.globals.maxSize = 10000 * 1024^2)
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

# set parameters ----------------------------------------------------------
logfc.threshold.run <- 0.25
min.pct.run <- 0.1
min.diff.pct.run <- 0
mc_name = "MC8"
mc_name = "MC2"

# FindClusters -----------------------------------------------------------------
srat <- FindClusters(srat, resolution = 0.5)
cat("Finished FindClusters!\n")

# Preprocess ---------------------------------------------------------------
Idents(srat) = "orig.ident"
srat@meta.data[,"MC_name"] = paste0("MC", as.numeric(srat@meta.data[,"integrated_snn_res.0.5"]))
table(srat@meta.data$MC_name)

# Findmarkers -------------------------------------------------------------
markers_df <- FindMarkers(object = srat, group.by = "MC_name", ident.1 = mc_name,
                          test.use = "wilcox", only.pos = F,
                          min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, 
                          verbose = T)
markers_df$gene_symbol <- rownames(markers_df)
path_markers <- paste0(dir_out, mc_name, ".notonlyposmarkers.logfcthreshold.", logfc.threshold.run, ".minpct.", min.pct.run, ".mindiffpct.", min.diff.pct.run, ".tsv")
write.table(x = markers_df, file = path_markers, quote = F, sep = "\t", row.names = F)
print("Finish writing the FindAllMarkers output!\n")

