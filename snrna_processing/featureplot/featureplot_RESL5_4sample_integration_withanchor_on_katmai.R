# Yige Wu @WashU Apr 2020
## for making featureplot for RESL5_4sample_integration

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
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
source("./ccRCC_drug_analysis/plotting.R")
library(Seurat)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# set dependencies --------------------------------------------------------
## set integration id
id_integration <- "RESL5_4sample_integration.withanchor.20200417.v1"
## set the path to the rds file for integrated object
path_rds <- "./Resources/Analysis_Results/snrna_processing/clustering/cluster_RESL5_4sample_integration_withanchor_on_katmai/20200417.v1/RESL5_4sample_integration.withanchor.20200416.v1.clustered.RDS"
## input RDS file
srat <- readRDS(file = path_rds)
## input marker gene table
gene2celltype_df <- fread("./Resources/Analysis_Results/dependencies/merge_celltypemarkergenes_btw_human_and_mouse/20200409.v1/celltypemarkergenes_mouse_human.rcc.20200409.v1.tsv", data.table = F)
## set the minimal % of cells expresssing the gene
min.exp.pct <- 0

# make plot ---------------------------------------------------------------
## make feature name
gene2celltype_df <- gene2celltype_df %>%
  mutate(feature_name = ifelse(Species == "Human", paste0("GRCh38-3.0.0.premrna-", Gene_Symbol), paste0("mm10-premrna---------", Gene_Symbol)))
## get feature names in RNA count data
featurenames <-  intersect(gene2celltype_df$feature_name, srat@assays$RNA@counts@Dimnames[[1]])
featurenames <- unique(featurenames)

for (featurename in featurenames) {
  p <- FeaturePlot(object = srat, features = featurename, 
                   min.cutoff = "q10", max.cutoff = "q90", sort.cell = TRUE,
                   cols = c("grey", "red"), reduction = "umap", label = T)
  
  # save output -------------------------------------------------------------
  file2write <- paste0(dir_out, featurename, ".RESL5_4sample_integration.withanchor.", run_id, ".png")
  png(filename = file2write, width = 3000, height = 2000, res = 150)
  print(p)
  dev.off()
}


