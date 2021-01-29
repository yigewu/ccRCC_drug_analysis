# Yige Wu @WashU Apr 2020
## plot dimplot with cell type annotated

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/plotting.R")
library(ggrastr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input cell type and umap data
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/snrna_processing/map_barcode/map_barcode2celltype_with_8sample_integration/20210129.v1/Barcode2CellType_UMAP.20210129.v1.tsv", data.table = F)

# make plot for all samples---------------------------------------------------------------
## make plot data
plotdata_df <- barcode2celltype_df %>%
  mutate(Cell_type_plot = Cell_Type.Short) %>%
  mutate(CellType_species = ifelse(Cell_type_plot != "Unknown", paste0(Cell_Type.Short, "-", Species_Cell), "Unknown"))
## make plot
p <- ggplot()
p <- p + geom_point_rast(data = subset(plotdata_df, CellType_species == "Unknown"), mapping = aes(x = UMAP_1, y = UMAP_2), size = 0.2, color = "grey90")
p <- p + geom_point_rast(data = subset(plotdata_df, CellType_species != "Unknown"), mapping = aes(x = UMAP_1, y = UMAP_2, color = CellType_species), size = 0.2, alpha = 0.8)
p <- p + theme_classic()
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + ggtitle(label = paste0("8 Sample Integration (All Cells)"), subtitle = paste0("CT & Cab & Sap & (Cab+Sap)"))
p
## wreite output
file2write <- paste0(dir_out, "umap_celltype_species.", ".png")
png(filename = file2write, width = 1200, height = 1000, res = 150)
print(p)
dev.off()

# make plot divided by each sample---------------------------------------------------------------
## make plot data
plotdata_df <- barcode2celltype_df %>%
  mutate(Cell_type_plot = Cell_Type.Short) %>%
  mutate(CellType_species = ifelse(Cell_type_plot != "Unknown", paste0(Cell_Type.Short, "-", Species_Cell), "Unknown")) %>%
  mutate(Treatment_group = str_split_fixed(string = Id_Sample, pattern = "-", n = 3)[,3]) %>%
  mutate(Treatment_group = gsub(x = Treatment_group, pattern = "2", replacement = "")) %>%
  mutate(Id_model = str_split_fixed(string = Id_Sample, pattern = "-", n = 3)[,1])
plotdata_df$Treatment_group <- factor(x = plotdata_df$Treatment_group, levels = c("CT", "Cabo", "Sap", "Cabo_Sap"))
## make plot
p <- ggplot()
p <- p + geom_point_rast(data = subset(plotdata_df, CellType_species == "Unknown"), mapping = aes(x = UMAP_1, y = UMAP_2), size = 0.2, color = "grey90")
p <- p + geom_point_rast(data = subset(plotdata_df, CellType_species != "Unknown"), mapping = aes(x = UMAP_1, y = UMAP_2, color = CellType_species), size = 0.2, alpha = 0.8)
p <- p + facet_grid(Id_model~Treatment_group)
p <- p + theme_bw(base_size = 15)
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + ggtitle(label = paste0("8 Sample Integration (All Cells)"))
p <- p + theme(legend.position = "bottom")
p
## wreite output
file2write <- paste0(dir_out, "umap_celltype_species.", "by_sample.",  ".png")
png(filename = file2write, width = 3000, height = 1500, res = 150)
print(p)
dev.off()
