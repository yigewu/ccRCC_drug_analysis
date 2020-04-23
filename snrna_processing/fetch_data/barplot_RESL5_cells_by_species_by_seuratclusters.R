# Yige Wu @WashU Apr 2020

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
id_integration <- "RESL5_4sample_integration.withanchor.20200417.v1"
## input fetched data
fetcheddata_df <- fread(input = "./Resources/Analysis_Results/snrna_processing/fetch_data/fetch_data_RESL5_4sample_integration_withanchor_on_katmai/20200423.v1/RESL5_4sample_integration.withanchor.20200417.v1.fetched_data.20200423.v1.tsv", data.table = F)

# make barplot ------------------------------------------------------------
## make plot data
plotdata_df <- fetcheddata_df %>%
  mutate(x = as.factor(seurat_clusters)) %>%
  mutate(Species = call)

## plot
p <- ggplot()
p <- p + geom_bar(data = plotdata_df, mapping = aes(x = x, fill = Species))
p <- p + xlab(label = "Id_Seurat_Cluster") + ylab(label = "Number_Cells")
p
# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "barplot.", id_integration, ".cells_by_species_by_seuratclusters.", run_id, ".png")
png(filename = file2write, width = 1000, height = 600, res = 150)
print(p)
dev.off()



