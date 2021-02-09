# Yige Wu @WashU Apr 2020
## for making dimplot for RESL5_4sample_integration

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
## input marker gene table
gene2celltype_df <- fread("./Resources/Knowledge/Gene_Lists/Kidney_Specific_EMT_Genes.20200911.v1.tsv", data.table = F)


# set parameters ----------------------------------------------------------
pct_thres <- 10
avgexp_thres <- 0

# get gene to plot --------------------------------------------------------
## make feature name
gene2celltype_df <- gene2celltype_df %>%
  mutate(feature_name = Gene)
## change ident
Idents(srat) <- "seurat_clusters"
## get feature names in RNA count data
featurenames <-  intersect(gene2celltype_df$feature_name, srat@assays$RNA@data@Dimnames[[1]])
featurenames <- unique(featurenames)
## get the pct expressed for each gene in each cluster
p <- DotPlot(object = srat, features = featurenames, col.min = 0, assay = "RNA")
expdata_df <- p$data
## transform the dataframe to matrix to better filter out genes with too low expressin
## filter genes based on the percentage expressed
pct_matrix <- dcast(data = expdata_df, formula = features.plot ~ id, value.var = "pct.exp")
genes_pct_filtered <- as.vector(pct_matrix[rowSums(pct_matrix[,unique(as.vector(expdata_df$id))] > pct_thres) >= 1, "features.plot"])
## filter genes based on the average expression
avgexp_matrix <- dcast(data = expdata_df, formula = features.plot ~ id, value.var = "avg.exp")
genes_exp_filtered <- as.vector(avgexp_matrix[rowSums(avgexp_matrix[,unique(as.vector(expdata_df$id))] > avgexp_thres) >= 1, "features.plot"])
## intersect
featurenames_filtered <- intersect(unique(genes_exp_filtered), unique(genes_pct_filtered))
print(length(featurenames_filtered))
cat("Finished making plot data!\n\n\n")

# save data ----------------------------------------------------------------
file2write <- paste0(dir_out, "Expression_Data.tsv")
write.table(x = expdata_df, file = file2write, sep = "\t", quote = F, row.names = F)

# plot not scaled -------------------------------------------------------------
cat("Start not scaled!\n\n\n")
plotdata_df <- expdata_df %>%
  filter(features.plot %in% featurenames_filtered)
expvalue_top <- quantile(x = plotdata_df$avg.exp, probs = 0.975)
plotdata_df <- plotdata_df %>%
  mutate(expvalue_plot = ifelse(avg.exp >= expvalue_top, expvalue_top, avg.exp))
plotdata_df$gene_group2 <- paste0(plyr::mapvalues(plotdata_df$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Gene_Group2), 
                                  "\n",
                                  "Markers")
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = features.plot, y = id, color = expvalue_plot, size = pct.exp), shape = 16)
p <- p + scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 9, name = "Spectral")[1:5]), guide = guide_legend(direction = "horizontal", nrow = 2, byrow = T))
p <- p + scale_size_continuous(range = c(0, 8), name="% Expressed", guide = guide_legend(direction = "horizontal"))
p <- p + facet_grid(.~gene_group2, scales = "free", space = "free", drop = T)
p <- p + theme(axis.text.x = element_text(angle = 90, face = "bold", size = 10))
p <- p + theme(axis.text.y = element_text( face = "bold", size = 12))
p <- p + theme(panel.spacing = unit(0, "lines"), panel.grid.major = element_line(colour = "grey80"), 
               panel.border = element_rect(color = "black", fill = NA, size = 0.5),
               panel.background = element_blank())
p <- p + theme(strip.background = element_rect(color = NA, fill = NA, size = 0.5), 
               strip.text.x = element_text(angle = 90, vjust = 0.5, size = 12, face = "bold"), strip.text.y = element_text(angle = 0, vjust = 0.5))
p <- p + theme(axis.title = element_blank())
p <- p + labs(colour = "Expression value")
p <- p + theme(legend.position = "bottom")
file2write <- paste0(dir_out, "EMTMarker.ExpNotScaled.png")
png(file = file2write, width = 900, height = 1500, res = 150)
print(p)
dev.off()
