# Yige Wu @WashU Apr 2020
## for making dimplot for RESL5_4sample_integration

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
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
library(Seurat)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# set dependencies --------------------------------------------------------
## set the path to the rds file for integrated object
path_rds <- "./Resources/Analysis_Results/snrna_processing/integration/run_humancells_8sample_integration_withanchor_on_katmai/20210208.v1/Humancells_8sample_integration.withanchor.20210208.v1.RDS"
## input RDS file
srat <- readRDS(file = path_rds)
## input manual cluster
barcode2cluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/map_barcode/map_barcode2idmanualcluster_humancells_integration_withanchor/20210211.v1/barcode2idmanualcluster_umapdata.20210211.v1.tsv")

# set parameters for findmarkers ------------------------------------------
logfc.threshold.run <- 0
min.pct.run <- 0.1
min.diff.pct.run <- 0.1
## spcify assay
assay_process <- "RNA"
cat(paste0("Assay: ", assay_process, "\n"))
cat("###########################################\n")

# add manual cluster ------------------------------------------------------
srat@meta.data$Id_Manual_Cluster <- mapvalues(x = rownames(srat@meta.data), from = barcode2cluster_df$Barcode_Integrated, to = as.vector(barcode2cluster_df$Id_Manual_Cluster))
cat(table(srat@meta.data$Id_Manual_Cluster))
Idents(srat) <- "Id_Manual_Cluster"

# run findallmarkers ------------------------------------------------------
markers_df <- FindAllMarkers(object = srat, test.use = "wilcox", only.pos = F,
                             min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T, assay = assay_process)
markers_df$deg_gene_symbol <- rownames(markers_df)
cat(paste0("Finished FindAllMarkers", "\n"))
cat("###########################################\n")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "FindAllMarkers.", "Wilcox.", "ByCluster.", run_id, ".tsv")
write.table(x = markers_df, file = file2write, sep = "\t", quote = F, row.names = F)
cat(paste0("Finished write.table", "\n"))
cat("###########################################\n")

