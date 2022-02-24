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
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
enricher_top_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/pathway/unite_ora_msigdb_H_CP_treated_vs_control_human_proteins/20220216.v1/ora_msigdb_H_CP.treated_vs_control.human.proteins.decreased.top.20220216.v1.tsv")

# make plot data -----------------------------------------------------------
plotdata_df <- enricher_top_df %>%
  mutate(pathway_name = str_split_fixed(string = ID, pattern = "_", n = 2)[,2]) %>%
  mutate(pathway_database = str_split_fixed(string = ID, pattern = "_", n = 2)[,1]) %>%
  mutate(y_plot = paste0(tolower(pathway_name), " (", tolower(pathway_database), ")")) %>%
  mutate(x_plot = generatio_num) %>%
  mutate(treatment = str_split_fixed(string = test, pattern = "_", n = 9)[,5]) %>%
  mutate(rank = ifelse(is.na(my_ranks), ">10", my_ranks))
plotdata_df$rank <- factor(x = plotdata_df$rank, levels = c(as.character(1:10), ">10"))
## order pathways
orderdata_df <- plotdata_df %>%
  filter(treatment == "CaboSap") %>%
  arrange(desc(x_plot))
plotdata_df$y_plot <- factor(x = plotdata_df$y_plot, levels = orderdata_df$y_plot)
color_red <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[1]
color_green <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[3]
color_yellow <- RColorBrewer::brewer.pal(n = 7, name = "Set2")[6]
colors_rank <- c(paste0("grey", 6*(0:9)), "white")
names(colors_rank) <- c(as.character(1:10), ">10")

# make plot ---------------------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot, fill = treatment, size = Count, color = rank, stroke = !is.na(my_ranks)), shape = 21)
p <- p + xlab("Gene ratio")
p <- p + scale_fill_manual(values = c("Cabo" = color_red, "Sap" = color_green, "CaboSap" = color_yellow))
p <- p + scale_color_manual(values = colors_rank)
p <- p + theme_bw()
p <- p + theme(axis.title.y = element_blank(), 
               axis.text.x = element_text(color = "black"),
               axis.text.y = element_text(color = "black"), legend.position = "right")
p <- p + guides(size = guide_legend(ncol = 2), color = guide_legend(ncol = 2))
p

file2write <- paste0(dir_out, "ora_msigdb_H_CP_pathways.png")
png(file2write, width = 1000, height = 800, res = 150)
print(p)
dev.off()
