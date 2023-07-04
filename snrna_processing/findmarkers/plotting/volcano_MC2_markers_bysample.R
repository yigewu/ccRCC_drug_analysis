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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_unique_markers_for_selected_clusters_katmai/20230618.v1/MC2.logfcthreshold.0.25.minpct.0.1.mindiffpct.0.tsv")
# senescence_df = readxl::read_xlsx(path = "~/Downloads/41467_2022_32552_MOESM4_ESM.xlsx")

# plot --------------------------------------------------------------------
sample_id = "RESL10F-12462-CT2"
for (sample_id in unique(deg_df$sample)) {
  plot_data_df = deg_df %>%
    filter(sample == sample_id) %>%
    mutate(log10p_val_adj = -log10(p_val_adj)) %>%
    mutate(log2FC = ifelse(avg_log2FC > 7.5, 7.5,
                           ifelse(avg_log2FC < -7.5, -7.5, avg_log2FC))) %>%
    mutate(Expression_change = ifelse(log10p_val_adj > -log10(0.05), 
                                      ifelse(avg_log2FC > log2(2), "up",
                                             ifelse(avg_log2FC < log2(0.5), "down", "ns")),
                                      "ns")) %>%
    mutate(text = ifelse(gene_symbol == "ICAM1", "ICAM1", NA))
  cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
  p <- ggplot(data = plot_data_df, mapping = aes(x = log2FC, y = log10p_val_adj))
  p <- p + geom_point(aes(color = Expression_change), alpha = 0.8, shape = 16) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed") + 
    geom_vline(xintercept = c(log2(0.5), log2(2)),
               linetype = "dashed")
  p <- p +  geom_label_repel(aes(label = text), min.segment.length = unit(0, 'lines'),
                             force = 2)
  p <- p +  scale_colour_manual(values = cols)
  p <- p + labs(x = "log2(fold change) for MC2 vs. other clusters",
                y = "-log10(adjusted p-value)",
                colour = "Expression\nchange")
  p <- p +  theme_bw() +  
    theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) 
  file2write <- paste0(dir_out, sample_id, ".pdf")
  pdf(file2write, width = 6, height = 5, useDingbats = F)
  print(p)
  dev.off()
}
