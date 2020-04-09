# Yige Wu @WashU Apr 2020
## make cell cycle gene table for human and mouse
## reference: https://github.com/hbctraining/scRNA-seq/blob/master/lessons/cell_cycle_scoring.md
## biomart reference: https://www.bioconductor.org/packages/devel/bioc/vignettes/biomaRt/inst/doc/biomaRt.html

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/load_pkgs.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/functions.R")
packages = c(
  "RCurl"
)

for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  library(package = pkg_name_tmp, character.only = T)
}
packages = c(
  'biomaRt'
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
  }
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## retreive human cell cycle genes
cc_human_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv") 
cell_cycle_human_genes <- read.csv(text = cc_human_file)
## retreive human cell cycle genes
cc_mouse_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv") 
cell_cycle_mouse_genes <- read.csv(text = cc_mouse_file)
## input the feature file
feature_names <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/Resources/snRNA_Processed_Data/Cell_Ranger/outputs/GRCh38-3.0.0.premrna_and_mm10_premrna/RESL10F-12462-CT2/outs/raw_feature_bc_matrix/features.tsv.gz", 
                       data.table = F, col.names = c("feature_name_ensembl_gene_id", "feature_name_gene_name", "assay_type"))

# query human genes -------------------------------------------------------
feature_names_human <- feature_names[grepl(pattern = "GRCh38-3.0.0.premrna", feature_names$feature_name_ensembl_gene_id), c("feature_name_ensembl_gene_id", "feature_name_gene_name")]
feature_names_human$ensembl_gene_id <- str_split_fixed(string = feature_names_human$feature_name_ensembl_gene_id, pattern = "_", n = 2)[,2]
feature_names_human$gene_name <- str_split_fixed(string = feature_names_human$feature_name_gene_name, pattern = "_", n = 2)[,2]

## merge gene symbol info with cell cycle phase info
cell_cycle_human_genes <- merge(x = cell_cycle_human_genes, 
                                y = feature_names_human[, c("gene_name", "ensembl_gene_id")],
                                by.x = c("geneID"), by.y = c("ensembl_gene_id"), all.x = T)

# query mouse genes -------------------------------------------------------
feature_names_mouse <- feature_names[grepl(pattern = "mm10_premrna", feature_names$feature_name_ensembl_gene_id), c("feature_name_ensembl_gene_id", "feature_name_gene_name")]
feature_names_mouse$ensembl_gene_id <- str_split_fixed(string = feature_names_mouse$feature_name_ensembl_gene_id, pattern = "_________", n = 2)[,2]
feature_names_mouse$gene_name <- str_split_fixed(string = feature_names_mouse$feature_name_gene_name, pattern = "_________", n = 2)[,2]

## merge gene symbol info with cell cycle phase info
cell_cycle_mouse_genes <- merge(x = cell_cycle_mouse_genes, 
                                y = feature_names_mouse[, c("gene_name", "ensembl_gene_id")],
                                by.x = c("geneID"), by.y = c("ensembl_gene_id"), all.x = T)
# merge human and mouse ---------------------------------------------------
cell_cycle_genes <- rbind(cell_cycle_mouse_genes %>%
                            mutate(species = "mouse"),
                          cell_cycle_human_genes %>%
                            mutate(species = "human"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "cell_cycle_human_mouse_genes.", run_id, ".tsv")
write.table(x = cell_cycle_genes, file = file2write, quote = F, sep = "\t", row.names = F)

