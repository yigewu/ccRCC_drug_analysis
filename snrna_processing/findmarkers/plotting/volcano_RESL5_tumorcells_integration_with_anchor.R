# Yige Wu @WashU Jun 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
source("./ccRCC_drug_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_RESL5_tumorcells_integration_withanchor_on_katmai/20200617.v1/FindMarkers.Wilcox.RESL5.Tumor_cells.integration.withanchor.20200507.v1.Treated_vs_CT..tsv")

# plot by each treated sample vs control highlight selected genes ----------------------------------
for (sampleid_treated in unique(deg_df$sampleid_group1)) {
  ## filter for data to plot
  plot_data_df <- deg_df %>%
    filter(sampleid_group1 == sampleid_treated) %>%
    rename(deg_feature_name = deg_gene_symbol) %>%
    mutate(Species = ifelse(grepl(x = deg_feature_name, pattern = "GRCh38"), "Human", "Mouse")) %>%
    filter(Species == "Human") %>%
    mutate(deg_gene_symbol = str_split_fixed(string = deg_feature_name, pattern = "GRCh38-3.0.0.premrna-", n = 2)[,2]) %>%
    mutate(text_label = ifelse(p_val_adj < 0.05 & deg_gene_symbol %in% c(genes_cabo_targets, genes_pi3k_mtor), deg_gene_symbol, NA)) %>%
    mutate(color_point = ifelse(p_val_adj > 0.05, "FDR>=0.05",
                                ifelse(!is.na(text_label), "FDR<0.05 (highlight)", "FDR<0.05"))) %>%
    mutate(Log10p_val_adj = -log10(x = p_val_adj))
  
  ## plot
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df, mapping = aes(x = avg_logFC, y = Log10p_val_adj, color = color_point), alpha = 0.6)
  p <- p + scale_color_manual(values = c("FDR<0.05 (highlight)" = "red", "FDR<0.05" = "black", "FDR>=0.05" = "grey80"))
  p <- p + geom_text_repel(data = plot_data_df, mapping = aes(x = avg_logFC, y = Log10p_val_adj, label = text_label), color = "red")
  p <- p + theme_bw()
  p <- p + ggtitle(label = paste0(sampleid_treated, " vs Control Tumor cells snRNA"))
  p
  ## save output
  file2write <- paste0(dir_out, sampleid_treated, "_vs_CT.", "Selected_Genes", ".png")
  png(filename = file2write, width = 1000, height = 800, res = 150)
  print(p)
  dev.off()
}



# plot by each treated sample vs control top genes ----------------------------------
for (sampleid_treated in unique(deg_df$sampleid_group1)) {
  ## filter for data to plot
  plot_data_df <- deg_df %>%
    filter(sampleid_group1 == sampleid_treated) %>%
    rename(deg_feature_name = deg_gene_symbol) %>%
    mutate(Species = ifelse(grepl(x = deg_feature_name, pattern = "GRCh38"), "Human", "Mouse")) %>%
    filter(Species == "Human") %>%
    mutate(deg_gene_symbol = str_split_fixed(string = deg_feature_name, pattern = "GRCh38-3.0.0.premrna-", n = 2)[,2]) %>%
    mutate(text_label = ifelse(p_val_adj < 0.001 & abs(avg_logFC) > 5, deg_gene_symbol, NA)) %>%
    mutate(color_point = ifelse(p_val_adj > 0.05, "FDR>=0.05",
                                ifelse(!is.na(text_label), "FDR<0.05 (highlight)", "FDR<0.05"))) %>%
    mutate(Log10p_val_adj = -log10(x = p_val_adj))

  ## plot
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df, mapping = aes(x = avg_logFC, y = Log10p_val_adj, color = color_point), alpha = 0.6)
  p <- p + scale_color_manual(values = c("FDR<0.05 (highlight)" = "red", "FDR<0.05" = "black", "FDR>=0.05" = "grey80"))
  p <- p + geom_text_repel(data = plot_data_df, mapping = aes(x = avg_logFC, y = Log10p_val_adj, label = text_label), color = "red")
  p <- p + theme_bw()
  p <- p + ggtitle(label = paste0(sampleid_treated, " vs Control Tumor cells"))
  p
  ## save output
  file2write <- paste0(dir_out, sampleid_treated, "_vs_CT.", ".png")
  png(filename = file2write, width = 1000, height = 800, res = 150)
  print(p)
  dev.off()
}

