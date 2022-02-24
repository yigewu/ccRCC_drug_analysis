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
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input DEGs
deg_raw_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/test/cor_1month_controlproteins_vs_relativetumorvolumevscontrol_bytreatment/20220126.v1/spearman.1month.controlproteins_vs_relativetumorvolume.20220126.v1.tsv")
enricher_top_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/pathway/unite_ora_msigdb_H_CP_treated_relativetumorvolume_human_protein_markers_incontrol/20220203.v1/ora_msigdb_H_CP_human_protein_markers_in_control.top.20220203.v1.tsv")

# preprocess -------------------------------------------------------------
deg_df <- deg_raw_df %>%
  filter(abundance_type == "total_protein")
enricher_top_df <- enricher_top_df %>%
  filter(Count >=3) %>%
  mutate(treatment = str_split_fixed(string = test, pattern = "_", n = 9)[,5]) %>%
  mutate(treatment_matched = ifelse(treatment == "CaboSap", "Cabo+ Sap", treatment)) %>%
  mutate(assoc_direction = str_split_fixed(string = test, pattern = "_", n = 9)[,8])

# plot --------------------------------------------------------------------
treatment_tmp <- "Sap"
for (treatment_tmp in c("Cabo", "Sap", "Cabo+ Sap")) {
  ## get genes to highlight
  genestrings_tmp <- enricher_top_df$geneID[enricher_top_df$treatment_matched == treatment_tmp]
  genes_list_tmp <- str_split(string = genestrings_tmp, pattern = "\\/")
  genes_tmp <- unlist(genes_list_tmp)
  
  ## make plot data
  plot_data_df <- deg_df %>%
    filter(treatment == treatment_tmp) %>%
    mutate(text_gene = str_split_fixed(string = ID, pattern = "_", n = 2)[,1]) %>%
    mutate(label_plot = ifelse(text_gene %in% genes_tmp, text_gene, NA)) %>%
    mutate(y_plot = -log10(p.value)) %>%
    mutate(x_plot = rho) %>%
    mutate(direction = ifelse(x_plot > 0, "positive", "negative"))
  
  ## ggplot without jitter
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot), alpha = 0.6, size = 0.5)
  p <- p + geom_text_repel(data = plot_data_df,
                           mapping = aes(x = x_plot, y = y_plot, label = label_plot, color = direction), 
                           # color = "black", 
                           force = 2,
                           segment.size = 0.2, segment.alpha = 0.5,
                           max.overlaps = Inf)
  p <- p + scale_color_manual(values = c("positive" = "#E41A1C", "negative" = "#377EB8"))
  p <- p + xlab("Spearman's rho") + ylab("-Log10(P-value)")
  p <- p + ggtitle(label = paste0("Proteins (baseline) associated with ", treatment_tmp, "-treated tumor volume"))
  p <- p + theme_bw() + theme(legend.position = "none", title = element_text(size = 10))
  # p <- p + ggtitle(label = paste0(sampleid1, " vs ", sampleid2, "\nMouse Endothelial Cells snRNA Expression"))
  ## save output
  file2write <- paste0(dir_out, treatment_tmp, ".relativetumorvolume.spearman_corr.", "total_protein", ".png")
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


