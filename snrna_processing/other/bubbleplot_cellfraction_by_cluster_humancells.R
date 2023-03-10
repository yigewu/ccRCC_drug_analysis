# Yige Wu @WashU May 2020
## plot dimplot with cluster name

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
## input cell type and umap data
barcode2cluster_df <- fread(input = "./Resources/Analysis_Results/snrna_processing/fetch_data/fetch_data_humancells_8sample_integration_withanchor_on_katmai/20210209.v1/HumanCells_8sample.umap_data.20210209.v1.tsv", data.table = F)

# make plot data ----------------------------------------------------------
plot_data_df <- barcode2cluster_df %>%
  group_by(orig.ident, seurat_clusters) %>%
  summarise(count_bysample_bycluster = n())
count_bysample_df <- barcode2cluster_df %>%
  group_by(orig.ident) %>%
  summarise(count_bysample = n())  
plot_data_df <- merge(x = plot_data_df, y = count_bysample_df, by = c("orig.ident"), all.x = T)
plot_data_df <- plot_data_df %>%
  mutate(frac_bysample_bycluster = (count_bysample_bycluster/count_bysample)) %>%
  mutate(id_model = str_split_fixed(string = orig.ident, pattern = "-", n = 3)[,1]) %>%
  mutate(treatment_group = str_split_fixed(string = orig.ident, pattern = "-", n = 3)[,3]) %>%
  mutate(treatment_group = gsub(pattern = "2", replacement = "", x = treatment_group)) %>%
  mutate(cluster = factor(seurat_clusters))

# make bubbleplot ---------------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = plot_data_df, mapping = aes(x = cluster, y = orig.ident, size = frac_bysample_bycluster))
p <- p + facet_grid(rows = vars(id_model), scales = "free")
p <- p + theme_classic()
p