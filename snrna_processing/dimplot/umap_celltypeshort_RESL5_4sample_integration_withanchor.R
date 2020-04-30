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
id_model = "RESL5"
## set integration id
id_integration <- "RESL5_4sample_integration.withanchor.20200417.v1"
## input cell type and umap data
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/snrna_processing/map_barcode/map_barcode2celltype/20200424.v1/barcode2celltype_umapdata.20200424.v1.tsv", data.table = F)

# make plot ---------------------------------------------------------------
## make plot data
plotdata_df <- barcode2celltype_df %>%
  filter(Id_Model == "RESL5")
## make plot
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_Type.Short), size = 0.2)
p <- p + facet_grid(.~call)
p

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "umap_celltypeshort.", id_integration, ".", run_id, ".png")
png(filename = file2write, width = 3000, height = 1000, res = 150)
print(p)
dev.off()

