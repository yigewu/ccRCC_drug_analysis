# Yige Wu @WashU Apr 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the average expression
exp_wide_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/average_expression/unite_expression/unite_sct_data_humancells_mousecells_bycelltypeshorter/20210406.v1/SCT.data.AverageExpression.20210406.v1.tsv")

# identify genes to plot -------------------------------------------------
genes_plot <- c("PXDN", "PDXP", "GLUD2", "SUSD2")

# make plot data ----------------------------------------------------------
plotdata_wide_df <- exp_wide_df %>%
  filter(genesymbol_human %in% genes_plot)
plotdata_df <- melt(data = plotdata_wide_df)
summary(plotdata_df$value)
x_cap <- 4
plotdata_df$cell_group <- mapvalues(x = plotdata_df$variable, 
                                    from = c("Tumor_cells", "Endothelial.cells", "Myofibroblasts", "Macrophages", "Fibroblasts"), 
                                    to = c("Tumor cells", "Endothelial cells", "Myofibroblasts", "Macrophages", "Fibroblasts"))
plotdata_df <- plotdata_df %>%
  mutate(x_plot = ifelse(value > x_cap, x_cap, value))
plotdata_df$y_plot <- factor(x = plotdata_df$genesymbol_human, levels = rev(genes_plot))
# ## make colors
# colors_cellgroup <- c("#E7298A", "#1B9E77", "#7570B3","#000000", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00","#CC79A7", "#B2DF8A", "#FB9A99","grey50")
# names(colors_cellgroup) <- c("Tumor cells", "Normal epithelial cells", "Immune others", "B-cells", "CD4+ T-cells", "CD8+ T-cells", "Macrophages", "DC", "NK cells","Endothelial cells", "Myofibroblasts", "Fibroblasts","Unknown")
## make colors
colors_cellgroup <- RColorBrewer::brewer.pal(n = 5, name = "Set1")
names(colors_cellgroup) <- c("Tumor cells", "Endothelial cells", "Macrophages", "Fibroblasts", "Myofibroblasts")

# plot --------------------------------------------------------------------
p <- ggplot(data = plotdata_df, mapping = aes(x = y_plot, y = x_plot, fill = cell_group, color = "white"))
p <- p + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.6), alpha = 0.7)
p <- p + scale_fill_manual(values = colors_cellgroup)
p <- p + scale_color_manual(values = c("white" = NA))
p <- p + theme_classic(base_size = 12)
p <- p + coord_flip()
p <- p + ylab("Normalized expression")
p <- p + theme(panel.grid.major.y = element_line(size=.1, color="black" ))
p <- p + scale_y_reverse()
p <- p + scale_x_discrete(position = "top")
p <- p + theme(legend.position = "left")
p <- p + theme(axis.text.y = element_text(size = 12), axis.title.y = element_blank())
p <- p + theme(axis.text.x = element_text(size = 12), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "first")))
p <- p + ggtitle(label = paste0("ccRCC PDX sn Expression"))
file2write <- paste0(dir_out, "TreatmentMarkers", ".png")
png(file2write, width = 800, height = 500, res = 150)
print(p)
dev.off()

