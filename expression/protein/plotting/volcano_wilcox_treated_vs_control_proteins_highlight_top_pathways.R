# Yige Wu @WashU Jan 2022

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
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
# deg_raw_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/test/wilcox_paired_diff_protein_treated_vs_control/20220209.v1/Wilcox.Paired.Each_Treated_Group_vs_CT.20220209.v1.tsv")
deg_raw_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/pathway/map_human_proteins_treated_vs_control_entrezids/20220209.v1/Wilcox.Paired.Each_Treated_Group_vs_CT.20220209.v1.tsv")
enricher_top_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/pathway/unite_ora_msigdb_H_CP_treated_vs_control_human_proteins/20220216.v1/ora_msigdb_H_CP.treated_vs_control.human.proteins.20220216.v1.tsv")
## input correlation result
cor_result_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/test/cor_1month_proteinschanges_vs_relativetumorvolumevscontrol_bytreatment/20220217.v1/spearman.1month.controlproteins_vs_relativetumorvolume.20220217.v1.tsv")

# preprocess -------------------------------------------------------------
deg_df <- deg_raw_df
enricher_top_df <- enricher_top_df %>%
  filter(!is.na(my_ranks)) %>%
  mutate(treatment = str_split_fixed(string = test, pattern = "_", n = 9)[,5]) %>%
  mutate(treatment_matched = ifelse(treatment == "CaboSap", "Cabo+ Sap", treatment)) %>%
  mutate(assoc_direction = str_split_fixed(string = test, pattern = "_", n = 9)[,6])
cor_result_df <- cor_result_df %>%
  mutate(genesymbol = str_split_fixed(string = ID, pattern = "_", n = 2)[,1])

# others ------------------------------------------------------------------
treatment_tmp <- "Cabo+ Sap"
for (treatment_tmp in c("Cabo", "Sap", "Cabo+ Sap")) {
  ## get genes to highlight
  genestrings_tmp <- enricher_top_df$geneID[enricher_top_df$treatment_matched == treatment_tmp]
  genes_list_tmp <- str_split(string = genestrings_tmp, pattern = "\\/")
  genes_highlight_tmp <- unique(unlist(genes_list_tmp))
  # genes_cor_tmp <- cor_result_df$genesymbol[cor_result_df$treatment == treatment_tmp & cor_result_df$p.value < 0.1 & cor_result_df$abundance_type == "total_protein"]
  genes_cor_tmp <- cor_result_df$genesymbol[cor_result_df$treatment == treatment_tmp & cor_result_df$p.value < 0.05 & cor_result_df$abundance_type == "total_protein"]
  
  genes_highlight_tmp <- genes_highlight_tmp[genes_highlight_tmp %in% genes_cor_tmp]
  
  ## make plot data
  plot_data_df <- deg_df %>%
    filter(group1 == treatment_tmp) %>%
    mutate(text_gene = PG.Gene) %>%
    mutate(y_plot = -log10(pvalue)) %>%
    # mutate(x_plot = diff_estimate) %>%
    mutate(x_plot = ifelse(diff_estimate > 2.5, 2.5,
                           ifelse(diff_estimate < -2.5, -2.5, diff_estimate))) %>%
    mutate(direction = ifelse(x_plot > 0, "positive", "negative")) %>%
    unique()
  
  plot_data_df <- plot_data_df %>%
    mutate(label_plot = ifelse(text_gene %in% genes_highlight_tmp & pvalue < 0.05, text_gene, NA))
  
  
  ## ggplot without jitter
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot), alpha = 0.6, size = 0.5)
  p <- p + geom_text_repel(data = subset(plot_data_df, x_plot > 0),
                           mapping = aes(x = x_plot, y = y_plot, label = label_plot, color = direction),
                           force = 2, xlim = c(1, 2.75),
                           segment.size = 0.2, segment.alpha = 0.5,
                           max.overlaps = Inf)
  p <- p + geom_text_repel(data = subset(plot_data_df, x_plot < 0),
                           mapping = aes(x = x_plot, y = y_plot, label = label_plot, color = direction),
                           force = 3, xlim = c(-2.75, -1),
                           segment.size = 0.2, segment.alpha = 0.5,
                           max.overlaps = Inf)
  p <- p + xlim(-2.5, 2.5)
  p <- p + scale_color_manual(values = c("positive" = "#E41A1C", "negative" = "#377EB8"))
  p <- p + xlab("Difference in log2(Intensity) for drug- vs. vehicle-treated tumor") + ylab("-Log10(P-value)")
  p <- p + ggtitle(label = paste0("Protein changes between ", treatment_tmp, "- vs. vehicle-treated tumor"))
  p <- p + theme_bw() + theme(legend.position = "none", title = element_text(size = 10))
  # p <- p + ggtitle(label = paste0(sampleid1, " vs ", sampleid2, "\nMouse Endothelial Cells snRNA Expression"))
  ## save output
  file2write <- paste0(dir_out, treatment_tmp, ".png")
  png(filename = file2write, width = 800, height = 600, res = 150)
  print(p)
  dev.off()
}

# ## ggplot with jitter
# jitter_obj <- position_jitter(h=0.1, w=0.1)
# p <- ggplot()
# p <- p + geom_point(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot), alpha = 0.6, size = 0.5, position=jitter_obj)
# p <- p + geom_text_repel(data = plot_data_df,
#                          mapping = aes(x = x_plot, y = y_plot, label = label_plot), color = "black", force = 2, alpha = 0.8, 
#                          position = jitter_obj,
#                          max.overlaps = Inf)
# p <- p + theme_bw()
# # p <- p + ggtitle(label = paste0(sampleid1, " vs ", sampleid2, "\nMouse Endothelial Cells snRNA Expression"))
# p

# # plot for combination plot--------------------------------------------------------------------
# treatment_tmp <- "Cabo+ Sap"
# genestrings_tmp <- enricher_top_df$geneID[enricher_top_df$treatment_matched == treatment_tmp &
#                                             enricher_top_df$ID %in% c("REACTOME_CELLULAR_SENESCENCE", "KEGG_LYSOSOME", 
#                                                                       "HALLMARK_MYC_TARGETS_V1")]
# genes_list_tmp <- str_split(string = genestrings_tmp, pattern = "\\/")
# genes_tmp <- unlist(genes_list_tmp); genes_tmp
# 
# ## make plot data
# plot_data_df <- deg_df %>%
#   filter(group1 == treatment_tmp) %>%
#   mutate(text_gene = PG.Gene) %>%
#   mutate(y_plot = -log10(pvalue)) %>%
#   mutate(x_plot = ifelse(diff_estimate > 2.5, 2.5,
#                          ifelse(diff_estimate < -2.5, -2.5, diff_estimate))) %>%
#   mutate(direction = ifelse(x_plot > 0, "positive", "negative"))
# 
# cutoff_x_high <- quantile(x = plot_data_df$x_plot, 0.9)
# cutoff_x_low <- quantile(x = plot_data_df$x_plot, 0.1)
# 
# plot_data_df <- plot_data_df %>%
#   mutate(label_plot = ifelse(text_gene %in% genes_tmp,text_gene, NA))
#   # mutate(label_plot = ifelse(text_gene %in% genes_tmp,
#   #                            ifelse(x_plot >= cutoff_x_high | x_plot <= cutoff_x_low, text_gene, NA), NA))
# 
# 
# p <- ggplot()
# p <- p + geom_point(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot), alpha = 0.6, size = 0.5)
# p <- p + geom_text_repel(data = subset(plot_data_df, x_plot > 0),
#                          mapping = aes(x = x_plot, y = y_plot, label = label_plot, color = direction), 
#                          # color = "black", 
#                          force = 2, xlim = c(1, 2.75),
#                          segment.size = 0.2, segment.alpha = 0.5,
#                          max.overlaps = Inf)
# p <- p + geom_text_repel(data = subset(plot_data_df, x_plot < 0),
#                          mapping = aes(x = x_plot, y = y_plot, label = label_plot, color = direction), 
#                          # color = "black", 
#                          force = 2, xlim = c(-2.75, -1),
#                          segment.size = 0.2, segment.alpha = 0.5,
#                          max.overlaps = Inf)
# p <- p + scale_color_manual(values = c("positive" = "#E41A1C", "negative" = "#377EB8"))
# p <- p + xlab("Spearman's rho") + ylab("-Log10(P-value)")
# p <- p + xlim(-2.5, 2.5)
# p <- p + ggtitle(label = paste0("Protein changes between ", treatment_tmp, "-treated vs. control"))
# p <- p + theme_bw() + theme(legend.position = "none", title = element_text(size = 10))
# # p <- p + ggtitle(label = paste0(sampleid1, " vs ", sampleid2, "\nMouse Endothelial Cells snRNA Expression"))
# ## save output
# file2write <- paste0(dir_out, treatment_tmp, ".highlighted.png")
# png(filename = file2write, width = 800, height = 600, res = 150)
# print(p)
# dev.off()


