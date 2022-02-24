# Yige Wu @ WashU 2022 Jan

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
source("./ccRCC_drug_analysis/plotting.R")
## set run id
version_tmp <- 4
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
enricher_top_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/pathway/unite_ora_msigdb_H_CP_treated_relativetumorvolume_human_protein_markers_incontrol/20220203.v1/ora_msigdb_H_CP_human_protein_markers_in_control.top.20220203.v1.tsv")

# make plot data -----------------------------------------------------------
plotdata_df <- enricher_top_df %>%
  filter(Count >=3) %>%
  mutate(y_plot = str_split_fixed(string = ID, pattern = "_", n = 2)[,2]) %>%
  mutate(number_genes_tested = str_split_fixed(string = GeneRatio, pattern = "\\/", n = 2)[,2]) %>%
  mutate(number_genes_tested = as.numeric(number_genes_tested)) %>%
  mutate(x_plot = Count/number_genes_tested) %>%
  mutate(treatment = str_split_fixed(string = test, pattern = "_", n = 9)[,5]) %>%
  mutate(assoc_direction = str_split_fixed(string = test, pattern = "_", n = 9)[,8])
plotdata_df$treatment <- factor(x = plotdata_df$treatment, levels = c("Cabo", "Sap", "CaboSap"))

# make plot ---------------------------------------------------------------
p <- ggplot()
p <- p + geom_bar(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot, fill = assoc_direction), stat = "identity")
p <- p + xlab("Gene ratio")
p <- p + scale_fill_manual(values = c("pos" = "#E41A1C", "neg" = "#377EB8"))
p <- p + facet_grid(treatment+assoc_direction~., scales = "free_y", space = "free_y")
p <- p + theme_bw()
p <- p + theme(axis.title.y = element_blank(), 
               axis.text.x = element_text(color = "black"),
               axis.text.y = element_text(color = "black"), legend.position = "none")
file2write <- paste0(dir_out, "ora_msigdb_H_CP_pathways.png")
png(file2write, width = 800, height = 800, res = 150)
print(p)
dev.off()
