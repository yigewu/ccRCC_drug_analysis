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
## set the path to the rds file for integrated object
path_rds <- "./Resources/Analysis_Results/snrna_processing/integration/run_humancells_8sample_integration_withanchor_on_katmai/20210208.v1/Humancells_8sample_integration.withanchor.20210208.v1.RDS"
## input RDS file
srat <- readRDS(file = path_rds)
DefaultAssay(srat) <- "RNA"
## set the minimal % of cells expresssing the gene
min.exp.pct <- 0
genes_process <- c("PXDN", "PDXP", "GLUD2", "SUSD2")

# get gene to plot --------------------------------------------------------
## change ident
Idents(srat) <- "orig.ident"
## get feature names in RNA count data
featurenames <-  intersect(genes_process, srat@assays$RNA@data@Dimnames[[1]])
featurenames <- unique(featurenames)
## get the pct expressed for each gene in each cluster
p <- DotPlot(object = srat, features = featurenames, col.min = 0, assay = "RNA")
expdata_df <- p$data
# print(expdata_df[1:4, 1:4])
## transform the dataframe to matrix to better filter out genes with too low expressin
plot_matrix <- dcast(data = expdata_df, formula = features.plot ~ id, value.var = "pct.exp")
## filter for genes that are expressed in >XX% (min.exp.pct) of one cluster at least
## replot with the filtered genes plus malignant cell marker genes
featurenames_filtered <- as.vector(plot_matrix[rowSums(plot_matrix[,unique(as.vector(expdata_df$id))] > min.exp.pct) >= 1, "features.plot"])
print(length(featurenames_filtered))
cat("Finished making plot data!\n\n\n")

# write expression data ---------------------------------------------------
file2write <- paste0(dir_out, "Expression.tsv")
write.table(x = expdata_df, file = file2write, sep = "\t", row.names = F, quote = F)
cat("Finished write.table!\n\n\n")

# make scaled dotplot ------------------------------------------------------------
cat("Start plotting scaled!\n\n\n")
p <- DotPlot(object = srat, features = featurenames_filtered, col.min = 0, assay = "RNA")
# p <- p + RotatedAxis()
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey50"),
               strip.text.x = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(angle = 90, size = 15, face = "bold"),
               strip.placement = "outside")
cat("Finished plotting scaled!\n\n\n")
file2write <- paste0(dir_out, "dotplot.", "scaled.", "png")
png(filename = file2write, width = 1000, height = 1000, res = 150)
print(p)
dev.off()
cat("Finished writing png!\n\n\n")

# plot not scaled -------------------------------------------------------------
cat("Starting not-scaled plot!\n\n\n")
plotdata_df <- expdata_df %>%
  filter(features.plot %in% featurenames_filtered)
expvalue_top <- quantile(x = plotdata_df$avg.exp, probs = 0.95)
plotdata_df <- plotdata_df %>%
  mutate(expvalue_plot = ifelse(avg.exp >= expvalue_top, expvalue_top, avg.exp))
## add facet
## plot
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = features.plot, y = id, color = expvalue_plot, size = pct.exp), shape = 16)
p <- p + scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 9, name = "Spectral")[1:5]), guide = guide_legend(direction = "horizontal", nrow = 2, byrow = T))
p <- p + scale_size_continuous(range = c(0, 8), name="% Expressed", guide = guide_legend(direction = "horizontal"))
p <- p + theme(axis.text.x = element_text(angle = 90, size = 10))
p <- p + theme(axis.text.y = element_text(size = 12))
p <- p + theme(panel.spacing = unit(0, "lines"), panel.grid.major = element_line(colour = "grey80"), 
               panel.border = element_rect(color = "black", fill = NA, size = 0.5),
               panel.background = element_blank())
p <- p + theme(strip.background = element_blank(),
               strip.text.x = element_text(angle = 0, vjust = 0.5))
p <- p + theme(axis.title = element_blank())
p <- p + labs(colour = "Expression value")
p <- p + theme(legend.position = "bottom")
cat("Finished plotting not-scaled!\n\n\n")
file2write <- paste0(dir_out, "dotplot.", "not_scaled.", "png")
png(filename = file2write, width = 800, height = 1000, res = 150)
print(p)
dev.off()
cat("Finished writing png!\n\n\n")



