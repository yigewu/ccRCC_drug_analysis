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
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
library(Seurat)
## set output directory
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
res_process <- c(0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4)

# FindClusters -----------------------------------------------------------------
srat <- FindClusters(srat, resolution = res_process)
cat("Finished FindClusters!\n")

# fetch data --------------------------------------------------------------
metadata_df <- srat@meta.data
metadata_df$barcode <- rownames(metadata_df)
umapdata_df <- FetchData(object = srat, vars = c("orig.ident", "UMAP_1", "UMAP_2"))
umapdata_df$barcode <- rownames(umapdata_df)
## merge
barcode_info_df <- merge(x = umapdata_df %>%
                           select(orig.ident, barcode, UMAP_1, UMAP_2),
                         y = metadata_df %>%
                           select(orig.ident, barcode, integrated_snn_res.0.1, integrated_snn_res.0.2, integrated_snn_res.0.3, integrated_snn_res.0.4, integrated_snn_res.0.5,
                                  integrated_snn_res.1, integrated_snn_res.2),
                         by = c("orig.ident", "barcode"), all.x = T)
## write output
file2write <- paste0(dir_out, "HumanCells.8sample.Metadata.ByResolution.", run_id, ".tsv")
write.table(x = barcode_info_df, file = file2write, quote = F, sep = "\t", row.names = F)
print("Finish writing the barcode info output!\n")

# Findmarkers -------------------------------------------------------------
markers_all_df <- NULL
for (resolution_tmp in res_process) {
  ## output to directory not seperated by run id, in case some of the iterations won't finish in one day
  path_markers <- paste0(dir_out_parent, "res.", resolution_tmp, ".tumorcellsreclustered.markers.logfcthreshold.", logfc.threshold.run, ".minpct.", min.pct.run, ".mindiffpct.", min.diff.pct.run, ".tsv")
  if (file.exists(path_markers)) {
    markers <- fread(data.table = F, input = path_markers)
    cat(paste0("Markers for resolution ", resolution_tmp, "exists, reading!\n"))
  } else {
    cat(paste0("Markers for resolution ", resolution_tmp, "doesn't exist, running FindMarkers!\n"))
    
    Idents(srat) <- paste0("integrated_snn_res.", resolution_tmp)
    print(table(Idents(srat)))
    
    markers_df <- FindAllMarkers(object = srat, test.use = "wilcox", only.pos = T,
                              min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T)
    markers_df$gene_symbol <- rownames(markers_df)
    markers_df$resolution <- resolution_tmp
    write.table(x = markers_df, file = path_markers, quote = F, sep = "\t", row.names = F)
  }
  markers_all_df <- rbind(markers_all_df, markers_df)
}
file2write <- paste0(dir_out, "markers_by_resolution.logfcthreshold.", 
                     logfc.threshold.run, ".minpct.", min.pct.run, ".mindiffpct.", min.diff.pct.run,
                     ".tsv")
write.table(x = markers_all_df, file = file2write, quote = F, sep = "\t", row.names = F)
print("Finish writing the FindAllMarkers output!\n")

