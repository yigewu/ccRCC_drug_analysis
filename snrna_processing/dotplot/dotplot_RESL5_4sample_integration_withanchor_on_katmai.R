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
DefaultAssay(srat) <- "RNA"
## input marker gene table
gene2celltype_df <- fread("./Resources/Analysis_Results/dependencies/merge_celltypemarkergenes_btw_human_and_mouse/20200409.v1/celltypemarkergenes_mouse_human.rcc.20200409.v1.tsv", data.table = F)
## set the minimal % of cells expresssing the gene
min.exp.pct <- 0

# get gene to plot --------------------------------------------------------
## make feature name
gene2celltype_df <- gene2celltype_df %>%
  mutate(feature_name = ifelse(Species == "Human", paste0("GRCh38-3.0.0.premrna-", Gene_Symbol), paste0("mm10-premrna---------", Gene_Symbol)))
## get feature names in RNA count data
featurenames <-  intersect(gene2celltype_df$feature_name, srat@assays$RNA@data@Dimnames[[1]])
featurenames <- unique(featurenames)
## get the pct expressed for each gene in each cluster
p <- DotPlot(object = srat, features = featurenames, col.min = 0)
plot_data <- p$data
## transform the dataframe to matrix to better filter out genes with too low expressin
plot_matrix <- dcast(data = plot_data, formula = features.plot ~ id, value.var = "pct.exp")
## filter for genes that are expressed in >XX% (min.exp.pct) of one cluster at least
## replot with the filtered genes plus malignant cell marker genes
featurenames_filtered <- as.vector(plot_matrix[rowSums(plot_matrix[,unique(as.vector(plot_data$id))] > min.exp.pct) >= 1, "features.plot"])
print(featurenames_filtered)
# make Dimplot ------------------------------------------------------------
p <- DotPlot(object = srat, features = featurenames_filtered, col.min = 0, assay = "RNA", split.by = "call")
p$data$species <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$feature_name, to = gene2celltype_df$Species)
p$data$gene_cell_type_group <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$feature_name, to = gene2celltype_df$Cell_Type_Group)
p$data$gene_cell_type1 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$feature_name, to = gene2celltype_df$Cell_Type1)
p$data$gene_cell_type2 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$feature_name, to = gene2celltype_df$Cell_Type2)
p$data$gene_cell_type3 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$feature_name, to = gene2celltype_df$Cell_Type3)
p$data$gene_cell_type4 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$feature_name, to = gene2celltype_df$Cell_Type4)
# p <- p + RotatedAxis()
p <- p + facet_grid(.~species + gene_cell_type_group + gene_cell_type1 + gene_cell_type2 + gene_cell_type3 + gene_cell_type4, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey50"),
               strip.text.x = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(angle = 90, size = 15, face = "bold"),
               strip.placement = "outside")

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "dotplot.", "RESL5_4sample_integration.withanchor.", run_id, ".png")
png(filename = file2write, width = 3000, height = 2000, res = 150)
print(p)
dev.off()


