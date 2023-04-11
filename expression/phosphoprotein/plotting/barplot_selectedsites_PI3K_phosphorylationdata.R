# Yige Wu @ WashU 2023 Apr

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

library(ggplot2)

# input data --------------------------------------------------------------
exp_df = fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Data_Files/phosphosite_matrix-log2_ratios-MD_norm(1).tsv")
metadata_df = fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Data_Files/Proteomics_Meta_lodaing_merged_02222022.tsv")

# filter samples ----------------------------------------------------------
metadata_process_df = metadata_df %>%
  filter(grepl(x = `Original ModelID`, "RESL")) %>%
  filter(`Original ModelID` != "RESL9")

# set plotting parameters -------------------------------------------------
gene_phosphosite_list = list()
gene_phosphosite_list[[1]] = c("AKT2", "S474")
gene_phosphosite_list[[2]] = c("MTOR", "S2448")
gene_phosphosite_list[[3]] = c("MTOR", "S2448S2450")
gene_phosphosite_list[[4]] = c("MTOR", "S2448S2454")
gene_phosphosite_list[[5]] = c("EIF4EBP1", "S65")
gene_phosphosite_list[[6]] = c("MAPK3", "T202Y204")

for (i in 1:length(gene_phosphosite_list)) {
  gene_plot = gene_phosphosite_list[[i]][1]
  phosphosite_plot= gene_phosphosite_list[[i]][2]
  
  # filter data -------------------------------------------------------------
  data_wide_df = exp_df %>%
    mutate(phosphosite = str_split_fixed(string = Phosphosite.Index, pattern="_", n=7)[,7]) %>%
    filter(Gene == gene_plot) %>%
    filter(phosphosite == phosphosite_plot)
  phosphosite_indexs = data_wide_df$Phosphosite.Index
  data_wide_df = data_wide_df[, metadata_process_df$`Sample aliquot ID`]
  for (j in 1:nrow(data_wide_df)) {
    data_long_df = as.data.frame(t(data_wide_df[j,]))
    data_long_df$aliquot_id = rownames(data_long_df)
    colnames(data_long_df) = c("value", "aliquot_id")
    plot_data_df = merge(data_long_df, metadata_process_df, by.x = c("aliquot_id"), by.y = c("Sample aliquot ID"), all.x = T)
    plot_data_df$value = as.numeric(plot_data_df$value)
    plot_data_df$value = plot_data_df$value - min(plot_data_df$value) + 0.01
    plot_data_df$x_group = factor(plot_data_df$`PIK3CA mutation`, levels = c("H1047R", "D350G", "WT"))
    plot_data_df = plot_data_df %>%
      arrange(Treatment) %>%
      mutate(Treatment_group = ifelse(Treatment == "", "None/Vehicle", Treatment)) %>%
      filter(!is.na(value))
    if (nrow(plot_data_df) < 1) {
      next()
    }
    plot_data_df$x_plot = factor(plot_data_df$`Original ModelID`, levels = plot_data_df$`Original ModelID`)
    
    
    # plot --------------------------------------------------------------------
    p <- ggplot()
    p <- p + geom_bar(data = plot_data_df, mapping = aes(x = x_plot, y = value, fill = Treatment_group), stat = "identity")
    p <- p + scale_fill_manual(values = c("None/Vehicle" = "grey50", "Cabozantinib" = "red", "Sapanisertib" = "green"))
    p <- p + facet_grid(~x_group, space = "free_x", scales = "free_x")
    p <- p + ylab(paste0(gene_plot, " ", phosphosite_plot, "  phosphorylation")) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    png(filename = paste0(dir_out, gene_plot, "_", phosphosite_indexs[j], ".png"), res = 150, width = 650, height = 650)
    print(p)
    dev.off() 
  }
}
