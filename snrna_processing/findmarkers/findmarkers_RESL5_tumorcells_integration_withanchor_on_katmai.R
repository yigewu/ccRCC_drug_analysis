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
path_rds <- "./Resources/Analysis_Results/snrna_processing/integration/run_RESL5_4sample_tumorcells_integration_withanchor_on_katmai/20200507.v1/RESL5.Tumor_cells.integration.withanchor.20200507.v1.RDS"
## set integration id
id_integration <- "RESL5.Tumor_cells.integration.withanchor.20200507.v1"
## input RDS file
srat <- readRDS(file = path_rds)
## set CT as group2
sampleids <- unique(srat@meta.data$orig.ident)
sampleid_group2 <- sampleids[grepl(x = sampleids, pattern = "-CT")]
sampleids_group1 <- sampleids[!grepl(x = sampleids, pattern = "-CT")]

# loop for each drug-treated sample ------------------------------------------------------
markers_all_df <- NULL
for (sampleid_group1 in sampleids_group1) {
  ## make new metadata
  metadata_tmp <- srat@meta.data
  metadata_tmp$integrated_barcode <- rownames(metadata_tmp)
  metadata_tmp <- metadata_tmp %>%
    mutate(group_findmarkers = ifelse(orig.ident == sampleid_group1, "group1",
                                      ifelse(orig.ident == sampleid_group2, "group2", "other")))
  srat@meta.data <- metadata_tmp
  rownames(srat@meta.data) <- metadata_tmp$integrated_barcode
  Idents(srat) <- "group_findmarkers"
  
  # count cells -------------------------------------------------------------
  cellcount_group_df <- metadata_tmp %>%
    select(group_findmarkers) %>%
    table() %>%
    as.data.frame() %>%
    rename(group_findmarkers = '.')
  
  # run findallmarkers ------------------------------------------------------
  markers_df <- FindMarkers(object = srat, test.use = "wilcox", ident.1 = "group1", ident.2 = "group2", logfc.threshold = 0)
  markers_df$deg_gene_symbol <- rownames(markers_df)
  markers_df$cellcount_group1_findmarkers <- cellcount_group_df$Freq[cellcount_group_df$group_findmarkers == "group1"]
  markers_df$cellcount_group2_findmarkers <- cellcount_group_df$Freq[cellcount_group_df$group_findmarkers == "group2"]
  markers_df$sampleid_group1 <- sampleid_group1
  markers_df$sampleid_group2 <- sampleid_group2
  markers_all_df <- rbind(markers_all_df, markers_df)
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "FindMarkers.", "Wilcox.", id_integration, ".Treated_vs_CT.", ".tsv")
write.table(x = markers_all_df, file = file2write, sep = "\t", quote = F, row.names = F)


