# Yige Wu @ WashU 2021 Mar

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
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
## input the protein data
exp_df <- fread("./Data_Freeze/v1.dataFreeze.washU_rcc/3.geneExp/v3.20210116/datafreeze.v3.kallisto.geneExp.protein_coding.tsv", data.table = F)
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/0.sample_info/v3.20210116/RCC_PDX_Samples.20210115.v2.xlsx")

# set parameters ----------------------------------------------------------
genes_filter <- c("ICAM1")
colnames_id <- colnames(exp_df)[!(colnames(exp_df) %in% sampleinfo_df$Analysis_ID)]

# plot --------------------------------------------------------------------
gene_plot <- genes_filter[1]
for (gene_plot in genes_filter) {
  ## make plot data
  plot_data_wide_df <- exp_df %>%
    filter(Name == gene_plot)
  plot_data_long_df <- melt(data = plot_data_wide_df, id.vars = colnames_id)
  plot_data_long_df$Treatment_length <- mapvalues(x = plot_data_long_df$variable, from = sampleinfo_df$Analysis_ID, to = as.vector(sampleinfo_df$Treatment.Month))
  plot_data_long_df$Treatment_length = paste0(plot_data_long_df$Treatment_length, " month")
  plot_data_long_df$Treatment.rnaformat <- mapvalues(x = plot_data_long_df$variable, from = sampleinfo_df$Analysis_ID, to = as.vector(sampleinfo_df$ShortTag))
  plot_data_long_df$Model_id <- mapvalues(x = plot_data_long_df$variable, from = sampleinfo_df$Analysis_ID, to = as.vector(sampleinfo_df$ModelID))
  plot_data_long_df <- as.data.frame(plot_data_long_df)
  plot_data_long_df$y_plot <- (plot_data_long_df$value - min(plot_data_long_df$value))/(max(plot_data_long_df$value) - min(plot_data_long_df$value))
  plot_data_long_df <- plot_data_long_df %>%
    dplyr::filter(Treatment.rnaformat %in% c("Control", "Treated.Cabo+Sap")) %>%
    dplyr::filter(Treatment_length == "1 month") %>%
    dplyr::mutate(Treatment = ifelse(Treatment.rnaformat == "Treated.Cabo+Sap", "Cabo+Sap", "Control"))
  plot_data_long_df$Treatment <- factor(x = plot_data_long_df$Treatment, levels = c("Control", "Cabo+Sap"))
  plot_data_long_df$Model_id <- factor(x =   plot_data_long_df$Model_id, levels = c("RESL5", "RESL10", "RESL11", "RESL3", "RESL12", "RESL4"))
  
  stat.test <- plot_data_long_df %>%
    filter(Model_id %in% c("RESL5", "RESL10")) %>%
    group_by(Model_id) %>%
    t_test(y_plot ~ Treatment, paired = F, 
           comparisons = list(c("Control", "Cabo+Sap"))) %>%
    adjust_pvalue() %>%
    add_significance("p")
  # stat.test <- stat.test %>% add_xy_position(fun = "mean_se", x = "Treatment")
  stat.test <- stat.test %>% add_xy_position(fun = "mean", x = "Treatment")
  
  ## plot
  p <- ggbarplot(data = plot_data_long_df, 
                 x = "Treatment", y = "y_plot", fill = "Treatment",
                 # add = "mean_se", 
                 add = c("mean", "dotplot"),
                 position = position_dodge(),
                 facet.by = "Model_id", nrow = 1)
  p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01)
  p <- p + labs(y = paste0(gene_plot, " gene expression\n(normalized)"))
  p <- p + theme_classic(base_size = 15)
  p <- p + theme(axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(),
                 legend.position = "bottom")
  p
  ## write output
  file2write <- paste0(dir_out, gene_plot, ".pdf")
  pdf(file2write, width = 5.5, height = 3, useDingbats = F)
  print(p)
  dev.off()
  file2write <- paste0(dir_out, gene_plot, ".png")
  png(file2write, width = 1500, height = 800, res = 150)
  print(p)
  dev.off()
}
