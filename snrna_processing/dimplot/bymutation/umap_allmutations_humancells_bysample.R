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
umap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/fetch_data/fetch_data_individual_sample_humancells_on_katmai/20210224.v1/HumanCells.Individual_sample.umap_data.20210224.v1.tsv")

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
snRNA_mutation_df <- snRNA_mutation_df %>%
  mutate(barcode_individual = str_split_fixed(string = barcode, pattern = "\\-", n = 3)[,1])

for (sampleid_tmp in unique(umap_df$orig.ident)) {
  ## input barcodes with mapped varaint alleles and reference alleles
  var_cells_df <- snRNA_mutation_df %>%
    filter(aliquot == sampleid_tmp) %>%
    filter(Is_Silent == F) %>%
    filter(allele_type == "Var")
  
  ref_cells_df <- snRNA_mutation_df %>%
    filter(aliquot == sampleid_tmp) %>%
    filter(Is_Silent == F) %>%
    filter(gene_symbol %in% var_cells_df$gene_symbol) %>%
    filter(allele_type == "Ref")
  
  ## make data frame for plotting
  umap_mut_df <- umap_df %>%
    filter(orig.ident == sampleid_tmp) %>%
    mutate(barcode_individual = barcode)
  ## merge with variant read info
  umap_mut_df <- merge(umap_mut_df, var_cells_df, by = c("barcode_individual"), all.x = T)
  
  ### create read type, distinguish variant allele and reference allele
  umap_mut_df$cell_allele_type <- "NA"
  umap_mut_df$cell_allele_type[!is.na(umap_mut_df$allele_type) & umap_mut_df$allele_type == "Var"] <- "Var"
  table(umap_mut_df$cell_allele_type)
  
  ### order the data frame so that cells mapped with variant allele will come on top
  umap_mut_df <- rbind(umap_mut_df[umap_mut_df$cell_allele_type == "NA",],
                       umap_mut_df[umap_mut_df$cell_allele_type == "Var",])
  
  ## add a column for drive genes
  umap_mut_df <- umap_mut_df %>%
    mutate(mutation_cat = ifelse(gene_symbol %in% ccRCC_drivers, "SMG",
                                 ifelse(gene_symbol == "C8orf76", "Shared", "Others"))) %>%
    mutate(Driver_Gene_Mutation = ifelse(gene_symbol %in% ccRCC_drivers, "TRUE", "FALSE")) %>%
    mutate(cell_allele_type_text = ifelse(cell_allele_type == "NA", "others", "cells with the variant read(s)"))
  ## make color palette for different read types
  colors_cell_allele_type <- c("#E31A1C", "grey70")
  names(colors_cell_allele_type) <- c("cells with the variant read(s)", "others")
  
  # plot all cells together -------------------------------------------------
  plot_data_df <- umap_mut_df
  p <- ggplot()
  p <- p + geom_point_rast(data = plot_data_df[plot_data_df$cell_allele_type == "NA",], mapping = aes(x = UMAP_1, y = UMAP_2), alpha = 0.8, size = 0.5, color = colors_cell_allele_type["others"])
  p <- p + geom_point_rast(data = plot_data_df[plot_data_df$cell_allele_type == "Var" & !(plot_data_df$gene_symbol %in% ccRCC_SMGs),], mapping = aes(x = UMAP_1, y = UMAP_2), alpha = 0.5, size = 0.8, color = colors_cell_allele_type["cells with the variant read(s)"])
  p <- p + geom_point_rast(data = plot_data_df[plot_data_df$cell_allele_type == "Var" & (plot_data_df$gene_symbol %in% ccRCC_SMGs),], mapping = aes(x = UMAP_1, y = UMAP_2), alpha = 0.8, size = 1, color = colors_cell_allele_type["cells with the variant read(s)"])
  p <- p + geom_text_repel(data = plot_data_df[!is.na(plot_data_df$gene_symbol),], 
                           mapping = aes(UMAP_1, UMAP_2, label = gene_symbol, colour = Driver_Gene_Mutation, size = Driver_Gene_Mutation))
  p <- p + scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey50"))
  p <- p + scale_size_manual(values = c("TRUE" = 4, "FALSE" = 3))
  p <- p +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  p <- p + theme(legend.position = "none")
  p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank())
  p
  file2write <- paste0(dir_out, sampleid_tmp, ".allmutations.png")
  png(file2write, width = 800, height = 900, res = 150)
  print(p)
  dev.off()
}
