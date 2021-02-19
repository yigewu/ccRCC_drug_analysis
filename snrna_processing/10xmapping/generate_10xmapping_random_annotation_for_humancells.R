# Yige Wu @WashU Feb 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the barcode info
umap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/fetch_data/fetch_data_humancells_8sample_integration_withanchor_on_katmai/20210209.v1/HumanCells_8sample.umap_data.20210209.v1.tsv")

# input seurat object -----------------------------------------------------
for (sampleid_tmp in unique(umap_df$orig.ident)) {
  anno_tab_tmp <- umap_df %>%
    filter(orig.ident == sampleid_tmp) %>%
    mutate(barcode_bam = gsub(pattern = "\\_.", replacement = "\\-1", x = barcode)) %>%
    select(barcode_bam) %>%
    mutate(random_group = "0")
  write.table(x = anno_tab_tmp, file = paste0(dir_out, sampleid_tmp, ".snRNA.", "AfterQC_Barcodes.tsv"), quote = F, row.names = F, sep = "\t", col.names = F)
}


