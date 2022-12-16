# Yige Wu @WashU Apr 2022
## plot cell type on integration UMAP

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
source("./ccRCC_drug_analysis//functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input DEGs
results_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_byres_humancells_8sample_integration_on_katmai/20220905.v1/markers_by_resolution.logfcthreshold.0.25.minpct.0.1.mindiffpct.0.tsv")

# make matrix data -----------------------------------------------------------------
plotdata_df <- results_df %>%
  filter(p_val_adj < 0.05) %>%
  mutate(diff_pct = (pct.1 - pct.2)) %>%
  filter(diff_pct >= 0.1) %>%
  # filter(diff_pct >= 0) %>%
  group_by(resolution, cluster) %>%
  summarise(number_degs = n())
plotdata_df$x_plot <- factor(plotdata_df$resolution)
textdata_df <- plotdata_df %>%
  group_by(resolution) %>%
  summarise(number_clusters = n()) %>%
  mutate(text = paste0("n=", number_clusters))
  # mutate(text = paste0("n=", number_clusters))
textdata_df$x_plot <- factor(textdata_df$resolution)


# plot --------------------------------------------------------------------
## set fontsize
fontsize_plot <- 24
pos <- position_jitter(width = 0.1, seed = 1)
y_cutoff <- 50
p <- ggplot()
p <- p + geom_boxplot(data = plotdata_df, mapping = aes(x = x_plot, y = number_degs, group = x_plot))
# p <- p + geom_point(data = plotdata_df, mapping = aes(x = x_plot, y = number_degs), position = pos, alpha = 0.8)
p <- p + geom_dotplot(data = plotdata_df, mapping = aes(x = x_plot, y = number_degs), 
                      color = NA, fill = "brown", binaxis='y', stackdir='center', dotsize=0.75, alpha = 0.75, shape = 16)
# p <- p + geom_hline(yintercept = 20, linetype = 2)
p <- p + geom_hline(yintercept = y_cutoff, linetype = 2)
p <- p + geom_text(data = textdata_df, mapping = aes(x = x_plot, y = 430, label = text), size = 8)
p <- p + scale_y_continuous(breaks = c(0, y_cutoff, 100, 200, 300, 400, 500))
p <- p + xlab("Resolution") + ylab("No. unique markers / cluster")
p <- p + theme_classic(base_size = fontsize_plot)
p <- p + theme(axis.text = element_text(color = "black", size = fontsize_plot), 
               axis.title = element_text(color = "black", size = fontsize_plot))
# save outputs ------------------------------------------------------------
file2write <- paste0(dir_out, "Number_of_DEGs_byresolution.pdf")
pdf(file2write, width = 11, height = 6, useDingbats = F)
print(p)
dev.off()

file2write <- paste0(dir_out, "Number_of_DEGs_byresolution.png")
png(file2write, width = 1000, height = 800, res = 150)
print(p)
dev.off()
