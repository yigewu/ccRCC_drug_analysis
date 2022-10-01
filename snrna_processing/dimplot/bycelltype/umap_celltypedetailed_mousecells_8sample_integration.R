# Yige Wu @WashU Apr 2020
## plot dimplot with cell type annotated

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2",
  "ggrastr"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
source("./ccRCC_drug_analysis//functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input cell type and umap data
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/snrna_processing/map_barcode/map_barcode2celltype_mousecells_integration_withanchor/20210225.v1/MouseCells.Barcode2CellType.20210225.v1.tsv", data.table = F)


# prepare -------------------------------------------------------------
# colors_celltype <- c(RColorBrewer::brewer.pal(n = 7, name = "Dark2"), RColorBrewer::brewer.pal(n = 3, name = "Set1")[2])
colors_celltype <- RColorBrewer::brewer.pal(n = 12, name = "Paired")[c(12, 6, 7, 5, 8, 1, 2, 10)]
names(colors_celltype) <- c("Endothelial cells", 
                            "Fibroblasts", "Fibroblasts Pdgfrb+", "Fibroblasts Col14a1+",  "Myofibroblasts",
                            "Macrophages", "Macrophages Tgfbi+", "Macrophages M2-like")
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(model = str_split_fixed(string = Id_Sample, pattern = "\\-", n = 3)[,1]) %>%
  mutate(model = gsub(x = model, pattern = '[A-Z]', replacement = "")) %>%
  mutate(model = paste0("RESL", model)) %>%
  mutate(treatment = str_split_fixed(string = Id_Sample, pattern = "\\-", n = 3)[,3]) %>%
  mutate(treatment = gsub(x = treatment, pattern = '[0-9]', replacement = ""))
barcode2celltype_df$treatment[barcode2celltype_df$treatment == "CT"] <- "control" 
barcode2celltype_df$treatment[barcode2celltype_df$treatment == "Sap"] <- "sapanisertib" 
barcode2celltype_df$treatment[barcode2celltype_df$treatment == "Cabo"] <- "cabozantinib" 
barcode2celltype_df$treatment[barcode2celltype_df$treatment == "Cabo_Sap"] <- "cabozantinib+sapanisertib" 
barcode2celltype_df$treatment <- factor(x = barcode2celltype_df$treatment, levels = c("control", "cabozantinib", "sapanisertib", "cabozantinib+sapanisertib"))


# make plot for all samples---------------------------------------------------------------
## make plot data
plotdata_df <- barcode2celltype_df %>%
  mutate(Cell_type_plot = Cell_Type.Detailed)
## make plot
p <- ggplot()
p <- p + geom_point_rast(data = subset(plotdata_df, Cell_type_plot == "Unknown"), mapping = aes(x = UMAP_1, y = UMAP_2), size = 0.2, color = "grey90")
p <- p + geom_point_rast(data = subset(plotdata_df, Cell_type_plot != "Unknown"), mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_type_plot), size = 0.2, alpha = 0.8)
p <- p + scale_color_manual(values = colors_celltype)
p <- p + guides(colour = guide_legend(override.aes = list(size=4), title = NULL, label.theme = element_text(size = 18), nrow = 4))
p <- p + theme_void()
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_blank())
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p <- p + theme(legend.position="bottom", aspect.ratio=1)
## write output
file2write <- paste0(dir_out, "umap_celltype.", "pdf")
pdf(file2write, width = 6, height = 6, useDingbats = F)
print(p)
dev.off()
# ## write output
# file2write <- paste0(dir_out, "umap_celltype.", "png")
# png(filename = file2write, width = 1200, height = 1000, res = 150)
# print(p)
# dev.off()

# make plot divided by each sample---------------------------------------------------------------
## make plot data
plotdata_df <- barcode2celltype_df %>%
  mutate(Cell_type_plot = Cell_Type.Detailed)
## make plot
p <- ggplot()
p <- p + geom_point_rast(data = subset(plotdata_df, Cell_type_plot == "Unknown"), mapping = aes(x = UMAP_1, y = UMAP_2), size = 0.2, color = "grey90")
p <- p + geom_point_rast(data = subset(plotdata_df, Cell_type_plot != "Unknown"), mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_type_plot), size = 0.2, alpha = 0.8)
p <- p + scale_color_manual(values = colors_celltype)
p <- p + facet_grid(rows = vars(model), cols = vars(treatment))
p <- p + guides(colour = guide_legend(override.aes = list(size=4), title = NULL, label.theme = element_text(size = 18), nrow = 3))
p <- p + theme_void(base_size = 15)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_blank())
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p <- p + theme(legend.position="none", aspect.ratio=1)
## write output
file2write <- paste0(dir_out, "umap_Cell_type_plot.", "by_sample.",  "pdf")
pdf(file2write, width = 10, height = 6, useDingbats = F)
print(p)
dev.off()
# ## write output
# file2write <- paste0(dir_out, "umap_Cell_type_plot.", "by_sample.",  ".png")
# png(filename = file2write, width = 3000, height = 1500, res = 150)
# print(p)
# dev.off()
