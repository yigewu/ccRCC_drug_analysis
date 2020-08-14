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
## input DEGs
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_RESL5_vs_RESL10_tumorcells_on_katmai/20200807.v1/FindMarkers.Wilcox.RESL.Tumor_cells.integration.withanchor.20200507.v1.RESL5_vs_RESL10..tsv")
## input genes to search for
genes_search_df <- fread(data.table = F, input = "../ccRCC_snRNA/Resources/Analysis_Results/dependencies/write_markers_by_intratumorheterogeneity_types/20200504.v1/markergenes_by_intratumorheterogeneity_types.20200504.v1.tsv")

# plot by each treated sample vs control highlight selected genes ----------------------------------
genes_search <- unique(genes_search_df$gene_symbol)

for (sampleid1 in unique(deg_df$sampleid_group1)) {
  sampleid2 <- unique(deg_df$sampleid_group2[deg_df$sampleid_group1 == sampleid1])
  ## filter for data to plot
  plot_data_df <- deg_df %>%
    dplyr::filter(sampleid_group1 == sampleid1) %>%
    rename(deg_feature_name = deg_gene_symbol) %>%
    mutate(Species = ifelse(grepl(x = deg_feature_name, pattern = "GRCh38"), "Human", "Mouse")) %>%
    filter(Species == "Human") %>%
    mutate(deg_gene_symbol = str_split_fixed(string = deg_feature_name, pattern = "GRCh38-3.0.0.premrna-", n = 2)[,2]) %>%
    mutate(text_label = ifelse(p_val_adj < 0.05 & deg_gene_symbol %in% genes_search, deg_gene_symbol, NA)) %>%
    mutate(color_point = ifelse(p_val_adj > 0.05, "FDR>=0.05",
                                ifelse(is.na(text_label), "FDR<0.05",
                                       ifelse(avg_logFC > 0, "FDR<0.05 (up)", "FDR<0.05 (down)")))) %>%
    mutate(Log10p_val_adj = -log10(x = p_val_adj))
  
  plot_data_df$color_point <- factor(plot_data_df$color_point, levels = c("FDR>=0.05", "FDR<0.05", "FDR<0.05 (up)", "FDR<0.05 (down)"))
  plot_data_df <- plot_data_df %>%
    arrange(color_point)
  
  ## plot
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df, mapping = aes(x = avg_logFC, y = Log10p_val_adj, color = color_point), alpha = 0.6)
  p <- p + scale_color_manual(values = c("FDR<0.05 (up)" = "red", "FDR<0.05 (down)" = "blue", "FDR<0.05" = "black", "FDR>=0.05" = "grey80"))
  p <- p + geom_text_repel(data = plot_data_df, mapping = aes(x = avg_logFC, y = Log10p_val_adj, label = text_label, color = color_point))
  p <- p + theme_bw()
  p <- p + ggtitle(label = paste0(sampleid1, " vs ", sampleid2, " snRNA"))
  p
  ## save output
  file2write <- paste0(dir_out, sampleid1, ".vs.", sampleid2, "Selected_Genes", ".png")
  png(filename = file2write, width = 1000, height = 800, res = 150)
  print(p)
  dev.off()
}

# plot by each treated sample vs control all significant genes ----------------------------------
genes_search <- unique(genes_search_df$gene_symbol)

for (sampleid1 in unique(deg_df$sampleid_group1)) {
  sampleid2 <- unique(deg_df$sampleid_group2[deg_df$sampleid_group1 == sampleid1])
  ## filter for data to plot
  plot_data_df <- deg_df %>%
    dplyr::filter(sampleid_group1 == sampleid1) %>%
    rename(deg_feature_name = deg_gene_symbol) %>%
    mutate(Species = ifelse(grepl(x = deg_feature_name, pattern = "GRCh38"), "Human", "Mouse")) %>%
    filter(Species == "Human") %>%
    mutate(deg_gene_symbol = str_split_fixed(string = deg_feature_name, pattern = "GRCh38-3.0.0.premrna-", n = 2)[,2]) %>%
    mutate(Log10p_val_adj = -log10(x = p_val_adj))
  y_showtext <- quantile(x = plot_data_df$Log10p_val_adj, probs = 0.9)
  
  plot_data_df <- plot_data_df %>%
    mutate(text_label = ifelse(p_val_adj < 0.05 & (deg_gene_symbol %in% genes_search | Log10p_val_adj >= y_showtext), deg_gene_symbol, NA)) %>%
    mutate(color_point = ifelse(p_val_adj > 0.05, "FDR>=0.05",
                                ifelse(is.na(text_label), "FDR<0.05",
                                       ifelse(avg_logFC > 0, "FDR<0.05 (up)", "FDR<0.05 (down)"))))
  
  plot_data_df$color_point <- factor(plot_data_df$color_point, levels = c("FDR>=0.05", "FDR<0.05", "FDR<0.05 (up)", "FDR<0.05 (down)"))
  plot_data_df <- plot_data_df %>%
    arrange(color_point)
  
  y_cap <- max(plot_data_df$Log10p_val_adj[!is.infinite(plot_data_df$Log10p_val_adj)])
  plot_data_df <- plot_data_df %>%
    mutate(y_capped = ifelse(Log10p_val_adj > y_cap, y_cap, Log10p_val_adj))
  y_limit <- y_cap*1.1
  
  ## plot
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df, mapping = aes(x = avg_logFC, y = y_capped, color = color_point), alpha = 0.6, size = 0.5)
  p <- p + scale_color_manual(values = c("FDR<0.05 (up)" = "red", "FDR<0.05 (down)" = "blue", "FDR<0.05" = "black", "FDR>=0.05" = "grey80"))
  p <- p + geom_text_repel(data = plot_data_df, mapping = aes(x = avg_logFC, y = y_capped, label = text_label, color = color_point), size = 3)
  p <- p + theme_bw()
  p <- p + ggtitle(label = paste0(sampleid1, " vs ", sampleid2, " snRNA"))
  p <- p + ylim(c(0, y_limit))
  p
  ## save output
  file2write <- paste0(dir_out, sampleid1, ".vs.", sampleid2,  ".png")
  png(filename = file2write, width = 1200, height = 800, res = 150)
  print(p)
  dev.off()
}
