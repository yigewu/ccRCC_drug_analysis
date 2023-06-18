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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/findmarkers/findmarkers_unique_markers_for_selected_clusters_katmai/20230413.v1/MC2.logfcthreshold.0.25.minpct.0.1.mindiffpct.0.tsv")
# senescence_df = readxl::read_xlsx(path = "~/Downloads/41467_2022_32552_MOESM4_ESM.xlsx")

# make plot data ----------------------------------------------------------
count_sig_df = deg_df %>%
  filter(p_val_adj < 0.05) %>%
  group_by(gene_symbol) %>%
  summarise(Num_sig = n())
plot_data_df = deg_df %>%
  mutate(log10p_val_adj = -log10(p_val_adj)) %>%
  group_by(gene_symbol) %>%
  summarise(log10p_val_adj = mean(log10p_val_adj),
            avg_log2FC = mean(avg_log2FC))
plot_data_df = merge(x = plot_data_df, y = count_sig_df, by = "gene_symbol", all.x = T)
plot_data_df$Num_sig[is.na(plot_data_df$Num_sig)] = 0
plot_data_df = plot_data_df %>%
  mutate(Expression_change = ifelse(log10p_val_adj > -log10(0.05) & Num_sig > 0, 
                                    ifelse(avg_log2FC > 0, "up", "down"),
                                    "ns"))

# plot --------------------------------------------------------------------
p <- ggplot(data = plot_data_df, mapping = aes(x = avg_log2FC, y = log10p_val_adj))
p <- p + geom_point(aes(color = Expression_change))
p
dev.off()
