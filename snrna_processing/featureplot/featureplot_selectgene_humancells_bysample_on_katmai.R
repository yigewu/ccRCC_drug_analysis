# Yige Wu @WashU Apr 2020
## for making featureplot for RESL5_4sample_integration

# set up libraries and output directory -----------------------------------
## set run id
version_tmp <- 3
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
library(Seurat)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# set dependencies --------------------------------------------------------
## set the path to the rds file for integrated object
path_rds <- "./Resources/Analysis_Results/snrna_processing/integration/run_humancells_8sample_integration_withanchor_on_katmai/20210208.v1/Humancells_8sample_integration.withanchor.20210208.v1.RDS"
## input DEGs
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/clusterprofiler/run_clusterprofiler_on_degs_humancells_bymanualcluster_filtered/20210216.v1/ORA.Result.20210216.v1.tsv")
## input RDS file
srat <- readRDS(file = path_rds)
DefaultAssay(srat) <- "RNA"

## input marker gene table
# genes_plot <- genes_rtk_cabo
# genes_plot <- c("VIM", "FN1", "MMP2", "CDH2", "HNF4G", "CA9")
# genes_plot <- c("EEF1A1", "EIF1", "RPS8", "RPS6", "HSPA1A", "HSPB1", "JUN", "DNAJB1",
#                 "MKI67", "RRM2", "BRCA1", "BRIP1", "ANLN",
#                 "MT-ND4", "MT-ND2", "MT-CYB", "MT-CO3",
#                 "IFIT2", "DDX58", "MX1", "IFI44L")
genes_plot <- sapply(deg_df$geneID[deg_df$DEG_Group == 3], function(x) {
  genes_tmp <- str_split(string = x, pattern = "\\/")[[1]]
  return(genes_tmp)
})
genes_plot <- unique(unlist(genes_plot))

## set the minimal % of cells expresssing the gene
min.exp.pct <- 0
# make plot ---------------------------------------------------------------
## get feature names in RNA count data
featurenames <-  intersect(genes_plot, srat@assays$RNA@counts@Dimnames[[1]])
featurenames <- unique(featurenames)

# modify srat object ------------------------------------------------------
print(dim(srat))

for (featurename in featurenames) {
  p <- FeaturePlot(object = srat, 
                   features = featurename,
                   split.by = "orig.ident",
                   min.cutoff = "q10", max.cutoff = "q90", sort.cell = TRUE,
                   cols = c("grey", "red"), reduction = "umap", label = T)
  
  # save output -------------------------------------------------------------
  file2write <- paste0(dir_out, featurename, ".featureplot", ".png")
  png(filename = file2write, width = 6000, height = 800, res = 150)
  print(p)
  dev.off()
}


