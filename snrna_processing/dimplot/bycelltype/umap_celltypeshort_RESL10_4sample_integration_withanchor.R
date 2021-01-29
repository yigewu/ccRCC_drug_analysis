# Yige Wu @WashU Apr 2020
## plot dimplot with cell type annotated

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
## set model id
id_model = "RESL10"
## set integration id
id_integration <- "RESL10_4sample_integration.withanchor.20200417.v1"
## input cell type and umap data
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/snrna_processing/map_barcode/map_barcode2celltype/20200617.v1/barcode2celltype_umapdata.20200617.v1.tsv", data.table = F)

# make plot for all samples---------------------------------------------------------------
## make plot data
plotdata_df <- barcode2celltype_df %>%
  filter(Id_Model == id_model) %>%
  mutate(CellType_Species = paste0(Cell_Type.Short, "-", Species_Cell))
## make plot
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = CellType_Species), size = 0.2, alpha = 0.8)
p <- p + theme_bw()
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + ggtitle(label = paste0(id_model, " Sample Integration (All Cells)"), subtitle = paste0("CT & Cab & Sap & (Cab+Sap)"))
p
## wreite output
file2write <- paste0(dir_out, "umap_celltype_species.", id_integration, ".", run_id, ".png")
png(filename = file2write, width = 1200, height = 1000, res = 150)
print(p)
dev.off()

# make plot divided by each sample---------------------------------------------------------------
## make plot data
plotdata_df <- barcode2celltype_df %>%
  filter(Id_Model == id_model) %>%
  mutate(CellType_Species = paste0(Cell_Type.Short, "-", Species_Cell))
unique(plotdata_df$Id_Sample)
plotdata_df$Id_Sample <- factor(x = plotdata_df$Id_Sample, levels = c("RESL10F-12462-CT2", "RESL10F-12465-Sap2", "RESL10F-12467-Cabo2", "RESL10F-12473-Cabo_Sap2"))
## make plot
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = CellType_Species), size = 0.2, alpha = 0.8)
p <- p + facet_grid(.~Id_Sample)
p <- p + theme_bw()
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + ggtitle(label = paste0(id_model, " Sample Integration (All Cells)"), subtitle = paste0("CT & Cab & Sap & (Cab+Sap)"))
p <- p + theme(legend.position = "bottom")
p
## wreite output
file2write <- paste0(dir_out, "umap_celltype_species.", id_integration, ".", "by_sample.", run_id, ".png")
png(filename = file2write, width = 3000, height = 1000, res = 150)
print(p)
dev.off()



