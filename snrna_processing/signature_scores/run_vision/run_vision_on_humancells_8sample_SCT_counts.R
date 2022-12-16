# Yige Wu @WashU Sep 2022
## https://yoseflab.github.io/VISION/articles/web_only/Seurat.html
## before run script, run conda activate vision

# set up libraries and output directory -----------------------------------
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
  "Seurat",
  "VISION"
)
if (!("VISION" %in% installed.packages()[,1])) {
  print(paste0(pkg_name_tmp, "is being installed!"))
  devtools::install_github("YosefLab/VISION", dependencies = T)
}
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_drug_analysis/functions.R")
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# dir_parent_out <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_Drug/Resources/Analysis_Results/snrna_processing/signature_scores/run_vision/run_vision_on_humancells_8sample_integrated/"
# dir.create(dir_parent_out)
# dir_out <- paste0(dir_parent_out, run_id, "/")
# dir.create(dir_out)

# input  ------------------------------------------------------
## input seurat object
path_rds <- "./Resources/Analysis_Results/snrna_processing/integration/run_humancells_8sample_integration_withanchor_on_katmai/20210208.v1/Humancells_8sample_integration.withanchor.20210208.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading the RDS file!\n")
## set paths to signature objects
signatures <- c("../ccRCC_snRNA/Resources/Knowledge/Databases/MSigDB/msigdb_v7.4_GMTs/h.all.v7.4.symbols.gmt",
                "../ccRCC_snRNA/Resources/Knowledge/Databases/MSigDB/msigdb_v7.4_GMTs/c2.cp.v7.4.symbols.gmt",
                "../ccRCC_snRNA/Resources/Knowledge/Databases/MSigDB/msigdb_v7.4_GMTs/c5.go.bp.v7.4.symbols.gmt")

# process -----------------------------------------------------------------
## this is important because the default assay for the seurat object is "integrated" and I think if not set to RNA assay Vision will just take integrated assay,
## which will give an error about  "rownames(data) = NULL. Expression matrix must have gene names as the rownames" because
### > rownames(srat$integrated@counts)
### NULL
### rownames(srat$integrated@data will give top variably expressed genes
# DefaultAssay(srat) <- "RNA"
options(mc.cores = 20)
## get counts data from the SCT assay
counts <- srat$SCT@counts
# Scale counts within each cell
n.umi <- colSums(counts)
scaled_counts <- t(t(counts) / n.umi) * median(n.umi)
## crate vision object
vision.obj <- Vision(scaled_counts, signatures = signatures, pool = F)
print("Finish creating the vision object!\n")
# Set the number of threads when running parallel computations
vision.obj <- analyze(vision.obj)
print("Finish analyze the vision object!\n")
sigScores <- getSignatureScores(vision.obj)
print("Finish getSignatureScores!\n")
file2write <- paste0(dir_out, "humancells.8sampleintegrated.Vision.sigScores.", run_id, ".RDS")
saveRDS(object = sigScores, file = file2write, compress = T)
sigCorr <- getSignatureAutocorrelation(vision.obj)
sigCorr$gene_set <- rownames(sigCorr)
print(head(sigCorr))
print("Finish getSignatureAutocorrelation!\n")
file2write <- paste0(dir_out, "humancells.8sampleintegrated.Vision.SignatureAutocorrelation.", run_id, ".tsv")
write.table(x = sigCorr, file = file2write, quote = F, sep = "\t", row.names = F)

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "humancells.8sampleintegrated.Vision.", run_id, ".RDS")
saveRDS(object = vision.obj, file = file2write, compress = T)


