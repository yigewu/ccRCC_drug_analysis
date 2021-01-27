# Yige Wu @ WashU 2021 Jan

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
source("./ccRCC_drug_analysis/plotting.R")
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the protein data
protein_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/preprocess/normalize_DIA_protein_data/20210111.v1/RCC_PDX.DIA_Protein.Log2.QuantileNormalized.20210111.v1.tsv")

# make data matrix --------------------------------------------------------
## Cabozantinib: VEGFR2 inhibitor, also inhibits c-Met, Ret, Kit, Flt-1/3/4, Tie2, and AXL
met_related_genes <- c("HGF", "MET", "AXL")
vegfr_genes <- c("FLT1", "KDR", "FLT3", "FLT4", "NRP1", "NRP2")
vegf_genes <- c("VEGFA", "VEGFB", "VEGFC", "VEGFD", "VEGFE")
other_cabo_related_genes <- c("KIT", "RET", "NTRK2", "TEK")

## sum up the genes
genes_plot <- c(met_related_genes, vegfr_genes, vegf_genes, other_cabo_related_genes)
genes_plot <- unique(genes2plot)
genes_plot <- c(genes_plot, "Flt1", "Flt4")
genes_plot[!(genes_plot %in% protein_df$PG.Genes)]
## filter by gene
plot_data_df <- protein_df[protein_df$PG.Genes %in% genes_plot,]
plot_data_mat <- as.matrix(plot_data_df[,5:ncol(plot_data_df)])
rownames(plot_data_mat) <- plot_data_df$PG.Genes

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color function for heatmap body colors
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
summary(as.vector(plot_data_mat))
colors_heatmapbody = circlize::colorRamp2(c(12, 
                                            16, 
                                            20), 
                                          c(color_blue, "white", color_red))

# plot  ------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = colors_heatmapbody,
                             na_col = color_na,
                             show_heatmap_legend = F)
p


# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "heatmap.png")
png(file2write, width = 1200, height = 600, res = 150)
print(p)
dev.off()
