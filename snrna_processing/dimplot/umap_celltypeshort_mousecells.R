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
## set integration id
id_integration <- "RESL.mouse_cells.integration.withanchor.20200814.v1"
## input the fetch data
fetchdata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/clustering/cluster_RESL_mousecells_integration_withanchor_on_katmai/20200814.v1/umap_data.RESL.mouse_cells.integration.withanchor.20200814.v1.PC40.Res0.5.tsv")
## input cell type and umap data
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/snrna_processing/map_barcode/map_barcode2celltype/20200617.v1/barcode2celltype_umapdata.20200617.v1.tsv", data.table = F)

# merge info --------------------------------------------------------------
merged_df <- merge(fetchdata_df %>%
                     rename(Id_Sample = orig.ident) %>%
                     mutate(Barcode = str_split_fixed(string = V1, pattern = "_", n = 2)[,1]), 
                   barcode2celltype_df %>%
                     mutate(Barcode = str_split_fixed(string = Barcode_Integrated, pattern = "_", n = 2)[,1]) %>%
                     dplyr::select(Id_Model, Id_Sample, Barcode, Cell_Type.Short, Species_Cell),
                   by.x = c("Id_Sample", "Barcode"),
                   by.y = c("Id_Sample", "Barcode"), all.x = T)
which(is.na(merged_df$Cell_Type.Short))

# make plot for all samples---------------------------------------------------------------
## make plot data
plotdata_df <- merged_df %>%
  mutate(CellType_Species = paste0(Cell_Type.Short, "-", Species_Cell))
## make plot
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = CellType_Species), size = 0.2, alpha = 0.8)
p <- p + theme_bw()
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + ggtitle(label = paste0("Mouse ", " Cell Integration"), subtitle = paste0("CT & Cab & Sap & (Cab+Sap)"))
p
## wreite output
file2write <- paste0(dir_out, "umap_celltype_species.", id_integration, ".", run_id, ".png")
png(filename = file2write, width = 1200, height = 1000, res = 150)
print(p)
dev.off()

# make plot divided by each sample---------------------------------------------------------------
## make plot data
plotdata_df <- merged_df %>%
  mutate(CellType_Species = paste0(Cell_Type.Short, "-", Species_Cell)) %>%
  mutate(Treatment = str_split_fixed(string = Id_Sample, pattern = "-", n = 3)[,3]) %>%
  mutate(Treatment = gsub(x = Treatment, pattern = '[0-9]', replacement = ""))
plotdata_df$Treatment <- factor(x = plotdata_df$Treatment, levels = c("CT", "Cabo", "Sap", "Cabo_Sap"))
plotdata_df$Id_Model <- factor(x = plotdata_df$Id_Model, levels = c("RESL5", "RESL10"))

## make plot
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = CellType_Species), size = 0.2, alpha = 0.8)
p <- p + facet_grid(Id_Model~Treatment)
p <- p + theme_bw()
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + ggtitle(label = paste0("Mouse ", " Cell Integration"), subtitle = paste0("CT & Cab & Sap & (Cab+Sap)"))
p <- p + theme(legend.position = "top")
p
## wreite output
file2write <- paste0(dir_out, "umap_celltype_species.", id_integration, ".", "by_sample.", run_id, ".png")
png(filename = file2write, width = 3000, height = 1500, res = 150)
print(p)
dev.off()
