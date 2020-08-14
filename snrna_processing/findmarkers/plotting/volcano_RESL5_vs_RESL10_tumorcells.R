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
genes_search_df <- readxl::read_excel(path = "./Resources/Knowledge/Gene_Lists/Targetable_Genes.20200812.xlsx")

# plot by each treated sample vs control highlight selected genes ----------------------------------
genes_search <- unique(genes_search_df$genesymbol)

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
    mutate(color_point = ifelse(p_val_adj > 0.05, "FDR>=0.05",
                                ifelse(avg_logFC > 0, "FDR<0.05 (up)", "FDR<0.05 (down)")))
  
  plot_data_df$color_point <- factor(plot_data_df$color_point, levels = c("FDR>=0.05", "FDR<0.05", "FDR<0.05 (up)", "FDR<0.05 (down)"))
  plot_data_df <- plot_data_df %>%
    arrange(color_point)
  
  ## find x value to show text
  x_high <- quantile(x = plot_data_df$avg_logFC[plot_data_df$avg_logFC > 0], prob = 0.90)
  x_high
  x_low <- quantile(x = plot_data_df$avg_logFC[plot_data_df$avg_logFC < 0], prob = 0.1)
  x_low
  plot_data_df <- plot_data_df %>%
    mutate(text_gene = ifelse((deg_gene_symbol %in% genes_search) & p_val_adj < 0.05, deg_gene_symbol, NA))
  
  ## cap y value
  y_cap <- max(plot_data_df$Log10p_val_adj[!is.infinite(plot_data_df$Log10p_val_adj)])
  plot_data_df <- plot_data_df %>%
    mutate(y_capped = ifelse(Log10p_val_adj > y_cap, y_cap, Log10p_val_adj))
  
  ## set y limit
  y_limit <- y_cap*1.1
  ## set x limit
  x_max <- max(plot_data_df$avg_logFC)
  x_limit_top <- x_max*1.1
  x_min <- min(plot_data_df$avg_logFC)
  x_limit_bottom <- x_min
  
  ## plot
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df, mapping = aes(x = avg_logFC, y = y_capped, color = color_point), alpha = 0.6, size = 0.5)
  p <- p + scale_color_manual(values = c("FDR<0.05 (up)" = "red", "FDR<0.05 (down)" = "blue", "FDR<0.05" = "black", "FDR>=0.05" = "grey80"))
  p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene) & avg_logFC > 0),
                           mapping = aes(x = avg_logFC, y = y_capped, label = text_gene), color = "black", force = 2, alpha = 0.8)
  p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene) & avg_logFC <= 0),
                           mapping = aes(x = avg_logFC, y = y_capped, label = text_gene), color = "black", force = 2)
  p <- p + theme_bw()
  p <- p + ggtitle(label = paste0(sampleid1, " vs ", sampleid2, " Tumor Cells snRNA Expression"))
  p <- p + ylim(c(0, y_limit))
  p <- p + xlim(c(x_limit_bottom, x_limit_top))
  
  p
  ## save output
  file2write <- paste0(dir_out, sampleid1, ".vs.", sampleid2, ".Targetable_Genes", ".png")
  png(filename = file2write, width = 1200, height = 800, res = 150)
  print(p)
  dev.off()
}


# plot by each treated sample vs control all significant genes ----------------------------------
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
    mutate(color_point = ifelse(p_val_adj > 0.05, "FDR>=0.05",
                                ifelse(avg_logFC > 0, "FDR<0.05 (up)", "FDR<0.05 (down)")))
  
  plot_data_df$color_point <- factor(plot_data_df$color_point, levels = c("FDR>=0.05", "FDR<0.05", "FDR<0.05 (up)", "FDR<0.05 (down)"))
  plot_data_df <- plot_data_df %>%
    arrange(color_point)
  
  ## find x value to show text
  x_high <- quantile(x = plot_data_df$avg_logFC[plot_data_df$avg_logFC > 0], prob = 0.90)
  x_high
  x_low <- quantile(x = plot_data_df$avg_logFC[plot_data_df$avg_logFC < 0], prob = 0.1)
  x_low
  plot_data_df <- plot_data_df %>%
    mutate(text_gene = ifelse((avg_logFC >= x_high | avg_logFC <= x_low) & p_val_adj < 0.05, deg_gene_symbol, NA))
  
  ## cap y value
  y_cap <- max(plot_data_df$Log10p_val_adj[!is.infinite(plot_data_df$Log10p_val_adj)])
  plot_data_df <- plot_data_df %>%
    mutate(y_capped = ifelse(Log10p_val_adj > y_cap, y_cap, Log10p_val_adj))
  
  ## set y limit
  y_limit <- y_cap*1.1
  ## set x limit
  x_max <- max(plot_data_df$avg_logFC)
  x_limit_top <- x_max*1.1
  x_min <- min(plot_data_df$avg_logFC)
  x_limit_bottom <- x_min
  
  ## plot
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df, mapping = aes(x = avg_logFC, y = y_capped, color = color_point), alpha = 0.6, size = 0.5)
  p <- p + scale_color_manual(values = c("FDR<0.05 (up)" = "red", "FDR<0.05 (down)" = "blue", "FDR<0.05" = "black", "FDR>=0.05" = "grey80"))
  p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene) & avg_logFC > 0),
                           mapping = aes(x = avg_logFC, y = y_capped, label = text_gene), color = "black", force = 2, alpha = 0.8)
  p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene) & avg_logFC <= 0),
                           mapping = aes(x = avg_logFC, y = y_capped, label = text_gene), color = "black", force = 2)
  p <- p + theme_bw()
  p <- p + ggtitle(label = paste0(sampleid1, " vs ", sampleid2, " Tumor Cells snRNA Expression"))
  p <- p + ylim(c(0, y_limit))
  p <- p + xlim(c(x_limit_bottom, x_limit_top))
  
  p
  ## save output
  file2write <- paste0(dir_out, sampleid1, ".vs.", sampleid2,  ".png")
  png(filename = file2write, width = 1200, height = 800, res = 150)
  print(p)
  dev.off()
}
