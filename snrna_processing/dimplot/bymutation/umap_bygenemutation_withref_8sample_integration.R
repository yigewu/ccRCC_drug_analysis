# Yige Wu @WashU Feb 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
source("./ccRCC_drug_analysis/plotting.R")
library(ggrastr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input 10xmapping result
snRNA_mutation_df <- fread("./Resources/Analysis_Results/snrna_processing/10xmapping/unite_10xmapping_humancells_heatmap_output/20210223.v1/10XMapping.HumanCells.mapping_heatmap.20210223.v1.tsv", data.table = F)
## input umap data
umap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/fetch_data/fetch_data_humancells_8sample_integration_withanchor_on_katmai/20210209.v1/HumanCells_8sample.umap_data.20210209.v1.tsv")

# check and remove the silent mutations ---------------------------------------------
snRNA_mutation_df <- snRNA_mutation_df %>%
  mutate(AA_Change = str_split_fixed(string = mutation, pattern = "-", n = 3)[,2]) %>%
  mutate(Ref_Allele = str_split_fixed(string = AA_Change, pattern = '[0-9]', n = 2)[,1])
snRNA_mutation_df$Alt_Allele <- sapply(X = snRNA_mutation_df$AA_Change, FUN = function(x) {
  text_vec <- str_split(string = x, pattern = '[0-9]')[[1]]
  return(text_vec[length(text_vec)])
})  
snRNA_mutation_df <- snRNA_mutation_df %>%
  mutate(Is_Silent = (Alt_Allele == Ref_Allele))
table(snRNA_mutation_df$Is_Silent)
# FALSE 
# 869

# preprocess --------------------------------------------------------------
## make color palette for different read types
colors_cell_allele_type <- c("#E31A1C", "#33A02C", "grey70")
names(colors_cell_allele_type) <- c("Var", "Ref", "NA")

## preprocess data
snRNA_mutation_df <- snRNA_mutation_df %>%
  mutate(barcode_individual = str_split_fixed(string = barcode, pattern = "\\-", n = 3)[,1])
var_cells_all_df <- snRNA_mutation_df %>%
  filter(Is_Silent == F) %>%
  filter(allele_type == "Var")
# for (gene_tmp in unique(var_cells_all_df$gene_symbol)) {
for (gene_tmp in "PIK3CA") {
    
  ## input barcodes with mapped varaint alleles and reference alleles
  ref_cells_df <- snRNA_mutation_df %>%
    filter(gene_symbol == gene_tmp) %>%
    filter(allele_type == "Ref")
  
  var_cells_df <- var_cells_all_df %>%
    filter(gene_symbol == gene_tmp)
  
  ## make data frame for plotting
  umap_mut_df <- umap_df %>%
    mutate(barcode_individual = str_split_fixed(string = barcode, pattern = "_", n = 3)[,1])
  ## merge with variant read info
  umap_mut_df <- merge(umap_mut_df, rbind(ref_cells_df, var_cells_df), by = c("barcode_individual"), all.x = T)
  
  ### create read type, distinguish variant allele and reference allele
  umap_mut_df <- umap_mut_df %>%
    mutate(cell_allele_type = ifelse(is.na(allele_type), "NA", allele_type))
  table(umap_mut_df$cell_allele_type)
  
  ### order the data frame so that cells mapped with variant allele will come on top
  umap_mut_df <- rbind(umap_mut_df[umap_mut_df$cell_allele_type == "NA",],
                       umap_mut_df[umap_mut_df$cell_allele_type == "Ref",],
                       umap_mut_df[umap_mut_df$cell_allele_type == "Var",])
  
  ## add a column for drive genes
  umap_mut_df <- umap_mut_df %>%
    mutate(cell_allele_type_text = cell_allele_type)
  
  # plot all cells together -------------------------------------------------
  plot_data_df <- umap_mut_df
  p <- ggplot()
  p <- p + geom_point_rast(data = plot_data_df[plot_data_df$cell_allele_type == "NA",], mapping = aes(x = UMAP_1, y = UMAP_2), alpha = 0.8, size = 0.5, color = colors_cell_allele_type["NA"])
  p <- p + geom_point_rast(data = plot_data_df[plot_data_df$cell_allele_type == "Ref",], mapping = aes(x = UMAP_1, y = UMAP_2), alpha = 0.8, size = 0.5, color = colors_cell_allele_type["Ref"])
  p <- p + geom_point_rast(data = plot_data_df[plot_data_df$cell_allele_type == "Var",], mapping = aes(x = UMAP_1, y = UMAP_2), alpha = 0.8, size = 1, color = colors_cell_allele_type["Var"])
  p <- p + geom_text_repel(data = plot_data_df[plot_data_df$cell_allele_type == "Var",], 
                           mapping = aes(UMAP_1, UMAP_2, label = gene_symbol), colour = "black", size = 4)
  p <- p +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  p <- p + theme(legend.position = "top")
  p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank())
  p
  file2write <- paste0(dir_out, gene_tmp, ".allsamples.mut.png")
  png(file2write, width = 800, height = 900, res = 150)
  print(p)
  dev.off()
}
