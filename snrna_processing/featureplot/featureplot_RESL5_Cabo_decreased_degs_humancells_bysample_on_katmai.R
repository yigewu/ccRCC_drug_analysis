# Yige Wu @WashU Apr 2020
## for making featureplot for RESL5_4sample_integration

# set up libraries and output directory -----------------------------------
## set run id
version_tmp <- "Tumor"
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
## input degs
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/examine_degs/filter_RESL5_Cab_decreased_degs_humancells/20210316.v1/DEGs_Decreased.Cabo_vs_CT.HumanCells.RESL5.ProteinFiltered.20210316.v1.tsv")
## input RDS file
srat <- readRDS(file = path_rds)
DefaultAssay(srat) <- "RNA"
srat$orig.ident <- factor(x = srat$orig.ident, levels = c("RESL10F-12462-CT2", "RESL10F-12465-Sap2", "RESL10F-12467-Cabo2", "RESL10F-12473-Cabo_Sap2",
                                                          "RESL5E-14541-CT2", "RESL5E-14542-Sap2", "RESL5E-14539-Cabo2", "RESL5E-14529-Cabo_Sap2"))

## input marker gene table
genes_plot <- deg_df$deg_gene_symbol

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


