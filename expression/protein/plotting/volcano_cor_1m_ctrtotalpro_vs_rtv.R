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
deg_raw_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/test/cor_1month_controlproteins_vs_relativetumorvolumevscontrol_bytreatment/20220126.v1/spearman.1month.controlproteins_vs_relativetumorvolume.20220126.v1.tsv")

# preprocess -------------------------------------------------------------
deg_df <- deg_raw_df %>%
  filter(abundance_type == "total_protein")

# plot --------------------------------------------------------------------
treatment_tmp <- "Sap"
for (treatment_tmp in c("Cabo", "Sap", "Cabo+ Sap")) {
  ## make plot data
  plot_data_df <- deg_df %>%
    filter(treatment == treatment_tmp) %>%
    mutate(text_gene = str_split_fixed(string = ID, pattern = "_", n = 2)[,1]) %>%
    mutate(label_plot = ifelse(p.value < 0.01, text_gene, NA)) %>%
    mutate(y_plot = -log10(p.value)) %>%
    mutate(x_plot = rho)
  ## ggplot without jitter
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot), alpha = 0.6, size = 0.5)
  p <- p + geom_text_repel(data = plot_data_df,
                           mapping = aes(x = x_plot, y = y_plot, label = label_plot), color = "black", force = 2, alpha = 0.8,
                           max.overlaps = Inf)
  p <- p + theme_bw()
  # p <- p + ggtitle(label = paste0(sampleid1, " vs ", sampleid2, "\nMouse Endothelial Cells snRNA Expression"))
  ## save output
  file2write <- paste0(dir_out, treatment_tmp, ".relativetumorvolume.spearman_corr.", "total_protein", ".png")
  png(filename = file2write, width = 800, height = 800, res = 150)
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


