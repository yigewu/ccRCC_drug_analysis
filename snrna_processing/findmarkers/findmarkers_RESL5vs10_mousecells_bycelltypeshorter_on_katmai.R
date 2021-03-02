# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## set run id
version_tmp <- 2
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
path_rds <- "./Resources/Analysis_Results/snrna_processing/integration/run_mousecells_8sample_integration_withanchor_on_katmai/20210208.v1/MouseCells_8sample_integration.withanchor.20210208.v1.RDS"
## input RDS file
srat <- readRDS(file = path_rds)
## input manual cluster group
barcode2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/map_barcode/map_barcode2celltype_mousecells_integration_withanchor/20210225.v1/MouseCells.Barcode2CellType.20210225.v1.tsv")

# set parameters ----------------------------------------------------------
## set ids
sampleids <- unique(srat@meta.data$orig.ident)
print(sampleids)
sampleids_group1 <- sampleids[grepl(x = sampleids, pattern = "RESL5")]
## spcify assay
assay_process <- "RNA"
DefaultAssay(srat) <- assay_process
cat(paste0("Assay: ", assay_process, "\n"))
cat("###########################################\n")
## set findmarkers parameters
logfc.threshold.run <- 0
min.pct.run <- 0.1
min.diff.pct.run <- 0.1

# loop for pair of RESL5 & RESL10 ------------------------------------------------------
markers_all_df <- NULL
for (sampleid_group1 in sampleids_group1) {
  print(paste0("Group1:", sampleid_group1))
  ## get RESL10 sample id
  treatment_name <- str_split_fixed(string = sampleid_group1, pattern = "-", n = 3)[,3]
  treatment_name <- gsub(x = treatment_name, pattern = '[0-9]', replacement = "")
  if (treatment_name == "Cabo_Sap") {
    sampleid_group2 <- sampleids[grepl(x = sampleids, pattern = "RESL10") & grepl(x = sampleids, pattern = treatment_name)]
  } else {
    sampleid_group2 <- sampleids[grepl(x = sampleids, pattern = "RESL10") & grepl(x = sampleids, pattern = treatment_name) & !grepl(x = sampleids, pattern = "Cabo_Sap")]
  }
  print(paste0("Group2:", sampleid_group2))
  ## make add manual cluster
  srat@meta.data$Id_Manual_Cluster <- mapvalues(x = rownames(srat@meta.data), from = barcode2celltype_df$Barcode_Integrated, to = as.vector(barcode2celltype_df$Cell_Type.Short))
  
  print(unique(srat@meta.data$Id_Manual_Cluster))
  for (cluster_tmp in unique(srat@meta.data$Id_Manual_Cluster)) {
    
    metadata_tmp <- srat@meta.data
    metadata_tmp$integrated_barcode <- rownames(metadata_tmp)
    metadata_tmp <- metadata_tmp %>%
      mutate(group_findmarkers = ifelse(Id_Manual_Cluster == cluster_tmp,
                                        ifelse(orig.ident == sampleid_group1, "group1",
                                               ifelse(orig.ident == sampleid_group2, "group2", "other")),
                                        "other"))
    rownames(metadata_tmp) <- metadata_tmp$integrated_barcode
    srat@meta.data <- metadata_tmp
    Idents(srat) <- "group_findmarkers"
    
    # count cells -------------------------------------------------------------
    cellcount_group_df <- metadata_tmp %>%
      select(group_findmarkers) %>%
      table() %>%
      as.data.frame() %>%
      rename(group_findmarkers = '.')
    
    # run findallmarkers ------------------------------------------------------
    markers_df <- FindMarkers(object = srat, test.use = "wilcox", 
                              min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T, assay = assay_process,
                              ident.1 = "group1", ident.2 = "group2")
    markers_df$deg_gene_symbol <- rownames(markers_df)
    markers_df$cellcount_group1_findmarkers <- cellcount_group_df$Freq[cellcount_group_df$group_findmarkers == "group1"]
    markers_df$cellcount_group2_findmarkers <- cellcount_group_df$Freq[cellcount_group_df$group_findmarkers == "group2"]
    markers_df$manual_cluster <- cluster_tmp
    markers_df$sampleid_group1 <- sampleid_group1
    markers_df$sampleid_group2 <- sampleid_group2
    markers_all_df <- rbind(markers_all_df, markers_df)
  }
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "FindMarkers.", "Wilcox.", ".RESL5_vs_RESL10.", "ByCelltypeshorter.", "tsv")
write.table(x = markers_all_df, file = file2write, sep = "\t", quote = F, row.names = F)


