# Yige Wu @WashU 2023

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_markers_for_selected_clusters_allmodels_katmai/20230618.v1/MC2.notonlyposmarkers.logfcthreshold.0.25.minpct.0.1.mindiffpct.0.tsv")
# senescence_df = readxl::read_xlsx(path = "~/Downloads/41467_2022_32552_MOESM4_ESM.xlsx")

# plot --------------------------------------------------------------------
plot_data_df = deg_df %>%
  mutate(log10p_val_adj = -log10(p_val_adj)) %>%
  mutate(plot_y = ifelse(log10p_val_adj > 350, 350, log10p_val_adj)) %>%
  mutate(plot_x = ifelse(avg_log2FC > 8, 8,
                         ifelse(avg_log2FC < -8, -8, avg_log2FC))) %>%
  mutate(Expression_change = ifelse(log10p_val_adj > -log10(0.05), 
                                    ifelse(avg_log2FC > log2(2), "up",
                                           ifelse(avg_log2FC < log2(0.5), "down", "ns")),
                                    "ns")) %>%
  mutate(text = ifelse(gene_symbol == "ICAM1", "ICAM1", NA))
cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
p <- ggplot(data = plot_data_df, mapping = aes(x = plot_x, y = plot_y))
p <- p + geom_point(aes(color = Expression_change), alpha = 0.8, shape = 16) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed")
p <- p +  geom_label_repel(aes(label = text), min.segment.length = unit(0, 'lines'),
                           force = 2, size = 6)
p <- p +  scale_colour_manual(values = cols)
p <- p + labs(x = "log2(fold change) for MC2 vs. other clusters",
              y = "-log10(adjusted p-value)",
              colour = "Expression change")
p <- p +  theme_bw(base_size = 20)
p <- p +  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        legend.position = "bottom") 
file2write <- paste0(dir_out, "MC2.markers.pdf")
pdf(file2write, width = 7, height = 5, useDingbats = F)
print(p)
dev.off()
