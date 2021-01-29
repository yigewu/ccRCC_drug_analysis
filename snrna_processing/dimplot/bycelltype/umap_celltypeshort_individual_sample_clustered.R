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
## input umap data by invidual sample clustered
umap_all_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/fetch_data/fetch_data_individual_sample_clustered_on_katmai/20210129.v1/Individual_sample_clustered.umap_data.20210129.v1.tsv")

# preprocess --------------------------------------------------------------
barcode2celltype_new_df <- merge(x = barcode2celltype_df %>%
                                   mutate(barcode = str_split_fixed(string = Barcode_Integrated, pattern = "_", n = "2")[,1]) %>%
                                   mutate(Cell_type_plot = Cell_Type.Short) %>%
                                   select(Id_Sample, barcode, Cell_type_plot, Species_Cell),
                                 y = umap_all_df,
                                 by.x = c("Id_Sample", "barcode"), by.y = c("orig.ident", "barcode"), all.x = T)


# make plot by samples---------------------------------------------------------------
for (id_sample_tmp in unique(barcode2celltype_df$Id_Sample)) {
  ## make plot data
  plotdata_df <- barcode2celltype_new_df %>%
    filter(Id_Sample == id_sample_tmp) %>%
    mutate(CellType_species = ifelse(Cell_type_plot != "Unknown", paste0(Cell_type_plot, "-", Species_Cell), "Unknown"))
  ## make plot
  p <- ggplot()
  p <- p + geom_point_rast(data = subset(plotdata_df, CellType_species == "Unknown"), mapping = aes(x = UMAP_1, y = UMAP_2), size = 0.2, color = "grey90")
  p <- p + geom_point_rast(data = subset(plotdata_df, CellType_species != "Unknown"), mapping = aes(x = UMAP_1, y = UMAP_2, color = CellType_species), size = 0.2, alpha = 0.8)
  p <- p + theme_classic()
  p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
  p <- p + ggtitle(label = paste0(id_sample_tmp, " (All Cells)"))
  p
  ## wreite output
  file2write <- paste0(dir_out, id_sample_tmp, ".umap_celltype_species.", ".png")
  png(filename = file2write, width = 1200, height = 1000, res = 150)
  print(p)
  dev.off()
}
