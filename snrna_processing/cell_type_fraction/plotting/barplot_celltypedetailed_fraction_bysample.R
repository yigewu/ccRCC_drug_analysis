# Yige Wu @WashU Mar 2021
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
source("./ccRCC_drug_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input fraction
celltype_frac_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/cell_type_fraction/calculate_celltypefraction_bycelltypedetailed_bysample/20220929.v1/CellTypeFraction.20220929.v1.tsv")

# prepare plot data -------------------------------------------------------
plot_data_df <-celltype_frac_df %>%
  mutate(Cell_Group_Plot = Cell_Type.Detailed) %>%
  mutate(model = str_split_fixed(string = Id_Sample, pattern = "\\-", n = 3)[,1]) %>%
  mutate(model = gsub(x = model, pattern = '[A-Z]', replacement = "")) %>%
  mutate(model = paste0("RESL", model)) %>%
  mutate(treatment = str_split_fixed(string = Id_Sample, pattern = "\\-", n = 3)[,3]) %>%
  mutate(treatment = gsub(x = treatment, pattern = '[0-9]', replacement = ""))
plot_data_df$treatment[plot_data_df$treatment == "CT"] <- "Control" 
plot_data_df$treatment[plot_data_df$treatment == "Sap"] <- "Sapanisertib" 
plot_data_df$treatment[plot_data_df$treatment == "Cabo"] <- "Cabozantinib" 
plot_data_df$treatment[plot_data_df$treatment == "Cabo_Sap"] <- "Cabozantinib+\nSapanisertib" 
plot_data_df$treatment <- factor(x = plot_data_df$treatment, levels = rev(c("Control", "Cabozantinib", "Sapanisertib", "Cabozantinib+\nSapanisertib")))
plot_data_df <- plot_data_df[order(plot_data_df$Id_Sample, by = plot_data_df$treatment), ]
plot_data_df$Id_Sample <- factor(x = plot_data_df$Id_Sample)
## make colors
colors_celltype <- c(RColorBrewer::brewer.pal(n = 12, name = "Paired")[c(12, 6, 7, 5, 8, 1, 2, 10)], RColorBrewer::brewer.pal(n = 4, name = "Dark2")[4])
names(colors_celltype) <- c("Endothelial cells", 
                            "Fibroblasts", "Fibroblasts Pdgfrb+", "Fibroblasts Col14a1+",  "Myofibroblasts",
                            "Macrophages", "Macrophages Tgfbi+", "Macrophages M2-like",
                            "Tumor cells")

# plot --------------------------------------------------------------------
textsize_plot <- 14
p <- ggplot()
p <- p + geom_col(data = plot_data_df, mapping = aes(x = treatment, y = Fraction_CellType_Sample, fill = Cell_Group_Plot), position = "stack")
p <- p + scale_fill_manual(values = colors_celltype)
p <- p + facet_grid(model~. , scales = "free", space = "free", drop = T)
p <- p + coord_flip()
p <- p + guides(fill = guide_legend(title = "Cell type", label.theme = element_text(size = textsize_plot), title.theme = element_text(size = textsize_plot)))
p <- p + ylab(label = "% cells by cell type")
p <- p + theme(legend.position = "right")
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), strip.background = element_blank(), 
               strip.text = element_text(color = "black", size = textsize_plot))
p <- p + theme(axis.title.y = element_blank(), 
               axis.title.x = element_text(size = textsize_plot), 
               axis.text = element_text(color = "black", size = textsize_plot), axis.ticks.y = element_blank())
file2write <- paste0(dir_out, "Cell_Group_Plot_Composition.", "pdf")
pdf(file2write, width = 6.75, height = 4, useDingbats = F)
print(p)
dev.off()
# file2write <- paste0(dir_out, "Cell_Group_Plot_Composition.", "png")
# png(file2write, width = 800, height = 400, res = 150)
# print(p)
# dev.off()

# plot --------------------------------------------------------------------
textsize_plot <- 14
for (celltype_tmp in unique(plot_data_df$Cell_Type.Detailed)) {
  p <- ggplot()
  p <- p + geom_col(data = subset(plot_data_df, Cell_Type.Detailed == celltype_tmp), mapping = aes(x = treatment, y = Fraction_CellType_Sample, fill = Cell_Group_Plot), position = "stack")
  p <- p + scale_fill_manual(values = colors_celltype)
  p <- p + facet_grid(model~. , scales = "free", space = "free", drop = T)
  p <- p + coord_flip()
  p <- p + guides(fill = guide_legend(title = "Cell type", label.theme = element_text(size = textsize_plot), title.theme = element_text(size = textsize_plot)))
  p <- p + ylab(label = "% cells by cell type")
  p <- p + theme(legend.position = "right")
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), strip.background = element_blank(), 
                 strip.text = element_text(color = "black", size = textsize_plot))
  p <- p + theme(axis.title.y = element_blank(), 
                 axis.title.x = element_text(size = textsize_plot), 
                 axis.text = element_text(color = "black", size = textsize_plot), axis.ticks.y = element_blank())
  file2write <- paste0(dir_out, celltype_tmp,  ".pdf")
  pdf(file2write, width = 6.75, height = 4, useDingbats = F)
  print(p)
  dev.off()
}
