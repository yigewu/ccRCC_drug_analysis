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
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2",
  "Seurat"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
library(Seurat)
## set output directory
source("./ccRCC_drug_analysis/functions.R")
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# set dependencies --------------------------------------------------------
##
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/snrna_processing/map_barcode/map_barcode2celltype_mousecells_integration_withanchor/20210225.v1/MouseCells.Barcode2CellType.20210225.v1.tsv", data.table = F)
## input marker gene table
gene2celltype_df <- fread("./Resources/Knowledge/Gene_Lists/Cell_Type_Marker_Genes/Mouse.Gene2CellType.20220928.tsv", data.table = F)
## input RDS file
path_rds <- "./Resources/Analysis_Results/snrna_processing/integration/run_mousecells_8sample_integration_withanchor_on_katmai/20210208.v1/MouseCells_8sample_integration.withanchor.20210208.v1.RDS"
srat <- readRDS(file = path_rds)
DefaultAssay(srat) <- "RNA"
## set the minimal % of cells expresssing the gene
min.exp.pct <- 10
textsize_plot <- 17

# process seurat object ---------------------------------------------------
srat@meta.data$cell_group <- mapvalues(x = rownames(srat@meta.data), 
                                       from = barcode2celltype_df$Barcode_Integrated, 
                                       to = barcode2celltype_df$Cell_Type.Detailed)
srat@meta.data$cell_group <- factor(x = srat@meta.data$cell_group, levels = c("Macrophages", "Macrophages Tgfbi+", "Macrophages M2-like",
                                                                              "Endothelial cells", 
                                                                              "Fibroblasts", "Fibroblasts Col14a1+", "Fibroblasts Pdgfrb+",
                                                                              "Myofibroblasts"))
table(srat@meta.data$cell_group)
Idents(srat) <- "cell_group"

# get gene to plot --------------------------------------------------------
# ## make feature name
# gene2celltype_df <- gene2celltype_df %>%
#   mutate(feature_name = Gene) %>%
#   filter(!(Cell_Type_Group %in% c("Neprhon_Epithelium", "Urothelium")))
# ## change ident
# srat@meta.data$id_by_cluster_species <- paste0(srat@meta.data$seurat_clusters, "_", srat@meta.data$call)
# Idents(srat) <- "id_by_cluster_species"
# ## get feature names in RNA count data
# featurenames <-  intersect(gene2celltype_df$feature_name, srat@assays$RNA@data@Dimnames[[1]])
# featurenames <- unique(featurenames)
# ## get the pct expressed for each gene in each cluster
# p <- DotPlot(object = srat, features = featurenames, col.min = 0, assay = "RNA")
# expdata_df <- p$data
# # print(expdata_df[1:4, 1:4])
# ## transform the dataframe to matrix to better filter out genes with too low expressin
# plot_matrix <- dcast(data = expdata_df, formula = features.plot ~ id, value.var = "pct.exp")
# ## filter for genes that are expressed in >XX% (min.exp.pct) of one cluster at least
# ## replot with the filtered genes plus malignant cell marker genes
# featurenames_filtered <- as.vector(plot_matrix[rowSums(plot_matrix[,unique(as.vector(expdata_df$id))] > min.exp.pct) >= 1, "features.plot"])
# print(length(featurenames_filtered))
# cat("Finished making plot data!\n\n\n")

# make scaled dotplot ------------------------------------------------------------
featurenames_filtered <- c("Ptprc", "Tgfbi", "Mrc1", "F13a1", "Flt1", "Pecam1", "Tek", "Col3a1", "Dcn", "Pdgfrb",  "Col14a1", "Acta2", "Mylk")
cat("Start plotting scaled!\n\n\n")
p <- DotPlot(object = srat, features = featurenames_filtered, col.min = 0, assay = "RNA")
p$data$features.plot <- factor(x = p$data$features.plot, levels = featurenames_filtered)
# p$data$gene_cell_type_group <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type_Group)
# p$data$gene_cell_type1 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type1)
# p$data$gene_cell_type2 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type2)
# p$data$gene_cell_type3 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type3)
# p$data$gene_cell_type4 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type4)
# p <- p + RotatedAxis()
# p <- p + facet_grid(.~gene_cell_type_group + gene_cell_type1 + gene_cell_type2 + gene_cell_type3 + gene_cell_type4, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey80"),
               axis.text.x = element_text(angle = 90, size = textsize_plot, face = "italic", color = "black", vjust = 0.5, hjust = 1), 
               axis.text.y = element_text(size = textsize_plot, color = "black"), 
               axis.title = element_blank(), legend.text = element_text(size = textsize_plot), legend.title = element_text(size = textsize_plot),
               strip.placement = "outside")
cat("Finished plotting scaled!\n\n\n")
file2write <- paste0(dir_out, "dotplot.", "scaled.", "pdf")
pdf(file2write, width = 8, height = 5, useDingbats = F)
print(p)
dev.off()
cat("Finished writing pdf!\n\n\n")
# file2write <- paste0(dir_out, "dotplot.", "scaled.", "png")
# png(filename = file2write, width = 3000, height = 1500, res = 150)
# print(p)
# dev.off()
# cat("Finished writing png!\n\n\n")

# plot not scaled -------------------------------------------------------------
cat("Starting not-scaled plot!\n\n\n")
p <- DotPlot(object = srat, features = featurenames_filtered, col.min = 0, assay = "RNA")
plotdata_df <- p$data
# plotdata_df <- expdata_df %>%
#   filter(features.plot %in% featurenames_filtered)
expvalue_top <- quantile(x = plotdata_df$avg.exp, probs = 0.95)
plotdata_df <- plotdata_df %>%
  mutate(expvalue_plot = ifelse(avg.exp >= expvalue_top, expvalue_top, avg.exp))
## add facet
# plotdata_df$gene_cell_type_group <- plyr::mapvalues(plotdata_df$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type_Group)
# plotdata_df$gene_cell_type1 <- plyr::mapvalues(plotdata_df$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type1)
# plotdata_df$gene_cell_type2 <- plyr::mapvalues(plotdata_df$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type2)
# plotdata_df$gene_cell_type3 <- plyr::mapvalues(plotdata_df$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type3)
# plotdata_df$gene_cell_type4 <- plyr::mapvalues(plotdata_df$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type4)
## plot
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = features.plot, y = id, color = expvalue_plot, size = pct.exp), shape = 16)
p <- p + scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 9, name = "Spectral")[1:5]), guide = guide_legend(direction = "horizontal", nrow = 1, byrow = T))
p <- p + scale_size_continuous(range = c(0, 8), name="% Expressed", guide = guide_legend(direction = "horizontal"))
# p <- p + facet_grid(.~gene_cell_type_group + gene_cell_type1 + gene_cell_type2 + gene_cell_type3 + gene_cell_type4, scales = "free", space = "free", drop = T)
p <- p + theme(axis.text.x = element_text(angle = 90, size = textsize_plot))
p <- p + theme(axis.text.y = element_text(size = textsize_plot))
p <- p + theme(panel.spacing = unit(0, "lines"), 
               panel.grid.major = element_line(colour = "grey80"), 
               panel.border = element_rect(color = "black", fill = NA, size = 0.5),
               panel.background = element_blank())
p <- p + theme(strip.background = element_blank(),
               strip.text.x = element_text(angle = 0, vjust = 0.5))
p <- p + theme(axis.title = element_blank())
p <- p + labs(colour = "Expression value")
p <- p + theme(legend.position = "bottom")
cat("Finished plotting not-scaled!\n\n\n")
# file2write <- paste0(dir_out, "dotplot.", "not_scaled.", "png")
# png(filename = file2write, width = 3000, height = 1700, res = 150)
# print(p)
# dev.off()
# cat("Finished writing png!\n\n\n")
file2write <- paste0(dir_out, "dotplot.", "notscaled.", "pdf")
pdf(file2write, width = 8, height = 5, useDingbats = F)
print(p)
dev.off()
cat("Finished writing pdf!\n\n\n")


