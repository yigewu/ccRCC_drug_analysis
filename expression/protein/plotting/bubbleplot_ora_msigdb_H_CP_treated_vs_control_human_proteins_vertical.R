# Yige Wu @ WashU 2022 Jan

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
packages = c(
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg_name_tmp)
  }
  library(package = pkg_name_tmp, character.only = T)
}

# input dependencies ------------------------------------------------------
# enricher_top_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/pathway/unite_ora_msigdb_H_CP_treated_vs_control_human_proteins/20220216.v1/ora_msigdb_H_CP.treated_vs_control.human.proteins.decreased.top.20220216.v1.tsv")
enricher_top_df1 <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/pathway/filter_ora_treated_vs_control_ttest_human_protein_msigdb_HCP/20220323.v1/ora_msigdb_H_CP.treated_vs_control.human.proteins.increased.top.20220323.v1.tsv")
enricher_top_df2 <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/pathway/filter_ora_treated_vs_control_ttest_human_protein_msigdb_HCP/20220323.v1/ora_msigdb_H_CP.treated_vs_control.human.proteins.decreased.top.20220323.v1.tsv")

# make plot data -----------------------------------------------------------
enricher_top_df <- rbind(enricher_top_df1, enricher_top_df2)
plotdata_df <- enricher_top_df %>%
  mutate(pathway_name = str_split_fixed(string = ID, pattern = "_", n = 2)[,2]) %>%
  mutate(pathway_database = str_split_fixed(string = ID, pattern = "_", n = 2)[,1]) %>%
  # mutate(y_plot = paste0(tolower(pathway_name), "|", tolower(pathway_database))) %>%
  mutate(y_plot = paste0(tolower(pathway_name), "(", tolower(pathway_database), ")")) %>%
  mutate(x_plot = generatio_num) %>%
  mutate(treatment = group1) %>%
  mutate(rank = ifelse(is.na(my_ranks), ">5", my_ranks))
  # mutate(rank = ifelse(is.na(my_ranks), ">10", my_ranks))
# plotdata_df$rank <- factor(x = plotdata_df$rank, levels = c(as.character(1:10), ">10"))
plotdata_df$rank <- factor(x = plotdata_df$rank, levels = c(as.character(1:5), ">5"))

## order pathways
orderdata_df <- plotdata_df %>%
  filter(treatment == "Cabo+ Sap") %>%
  # filter(treatment == "CaboSap") %>%
  arrange(desc(x_plot))
plotdata_df$y_plot <- factor(x = plotdata_df$y_plot, levels = orderdata_df$y_plot)
color_red <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[1]
color_green <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[3]
color_yellow <- RColorBrewer::brewer.pal(n = 7, name = "Set2")[6]
colors_rank <- c(paste0("grey", 6*(0:9)), "white")
names(colors_rank) <- c(as.character(1:10), ">10")
colors_rank <- c(paste0("grey", 6*(0:4)), "white")
names(colors_rank) <- c(as.character(1:5), ">5")

# make plot ---------------------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot, fill = treatment, size = Count, color = rank, stroke = !is.na(my_ranks)), shape = 21)
p <- p + facet_wrap(diff_direction~., scales = "free", nrow = 2)
p <- p + xlab("Gene ratio")
p <- p + scale_fill_manual(values = c("Cabo" = color_red, "Sap" = color_green, "Cabo+ Sap" = color_yellow))
p <- p + scale_color_manual(values = colors_rank)
p <- p + theme_bw()
p <- p + theme(axis.title.y = element_blank(), 
               axis.text.x = element_text(color = "black"),
               axis.text.y = element_text(color = "black"), legend.position = "right")
p <- p + guides(size = guide_legend(ncol = 2), color = guide_legend(ncol = 2))


# write output ------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_drug_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

file2write <- paste0(dir_out, "ora_msigdb_H_CP_pathways.pdf")
pdf(file2write, width = 7, height = 7, useDingbats = F)
print(p)
dev.off()
# file2write <- paste0(dir_out, "ora_msigdb_H_CP_pathways.png")
# png(file2write, width = 1000, height = 800, res = 150)
# print(p)
# dev.off()
