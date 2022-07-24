# Yige Wu @ WashU 2021 Nov

# set up libraries and output directory -----------------------------------
## set working directory
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
## input fold changes
foldchanges_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/get_foldchange_combo_vs_other_byprotein/20220317.v1/Protein_diff_combo_vs_other.20220317.v1.tsv")
## input the sample info
meta_data_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")

# prepare plotting parameters ---------------------------------------------
colors_treatment <- c("grey50", RColorBrewer::brewer.pal(n = 5, name = "Set1")[c(1,2,4,3)])
names(colors_treatment) <- c("Baseline", "Cabo", "Sap", "Cabo+ Sap", "Con")
# genes_process <- c("MFGE8", "PLA2G15", "SOD2", "CAVIN1", "LEMD3", "RPL27A", "MRPS6", "PTCD3")
genes_process <- c("PYCR1", "CLIC6")
genes_process <- c("PYCR1", "CHMP6", "MRM3", "ERO1B", "CLIC6", "COBLL1", "IGF2BP3", "DNAJC7")
cap_value <- 2
# cap_value <- 7

# plot all model -----------------------------------------------------------
treatment_tmp <- "Cabo"
treatment.month_tmp <- 1
model_ids_tmp <- unique(meta_data_control_df$ModelID)

## prepare plot data
plotdata_wide_df <- foldchanges_df[, c("PG.Gene", "PG.ProteinName", "PG.ProteinDescriptions", "PG.ProteinAccession", 
                                       paste0(model_ids_tmp, "_Combo.vs.", "Cabo", "_", treatment.month_tmp, "month"),
                                       paste0(model_ids_tmp, "_Combo.vs.", "Sap", "_", treatment.month_tmp, "month"),
                                       paste0(model_ids_tmp, "_Combo.vs.", "Con", "_", treatment.month_tmp, "month"))]
plotdata_wide_filtered_df <- plotdata_wide_df %>%
  filter(PG.Gene %in% genes_process)
# plotdata_wide_filtered_df <- plotdata_wide_filtered_df[rowSums(!is.na(plotdata_wide_filtered_df)) == ncol(plotdata_wide_filtered_df),]
plotdata_df <- plotdata_wide_filtered_df %>%
  melt() %>%
  filter(!is.infinite(value)) %>%
  mutate(comparison = str_split_fixed(string = variable, pattern = "_", n = 3)[,2]) %>%
  # filter(treatment == "Cabo+ Sap") %>%
  mutate(model_id = str_split_fixed(string = variable, pattern = "_", n = 3)[,1]) %>%
  mutate(log2FC = value) %>%
  mutate(x_plot = ifelse(log2FC > cap_value, cap_value,
                         ifelse(log2FC < (-cap_value), (-cap_value), log2FC)))
# orderdata_df <- plotdata_df %>%
#   filter(model_id == "RESL10") %>%
#   filter(treatment == "Sap") %>%
#   arrange(log2FC)
plotdata_df$model_id <- factor(x = plotdata_df$model_id, levels = rev(c("RESL11", "RESL10", "RESL5", "RESL12", "RESL3","RESL4")))
# plotdata_df$treatment <- factor(x = plotdata_df$treatment, levels = rev(c("Cabo", "Sap", "Cabo+ Sap")))
# plotdata_df$y_plot <- factor(x = plotdata_df$PG.Gene, levels = unique(orderdata_df$PG.Gene))
plotdata_df$y_plot <- factor(x = plotdata_df$PG.Gene, levels = rev(genes_process))

## plot
p <- ggplot(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot))
p <- p + geom_vline(xintercept = 0, linetype = 2)
p <- p + geom_bar(mapping = aes(fill = comparison), stat = "identity", position = 'dodge')
# p <- p + scale_fill_manual(values = colors_treatment)
p <- p + facet_grid(cols = vars(model_id))
p <- p + theme_classic()
p <- p + xlim(c(-cap_value, cap_value))
p <- p + xlab("Log2 (Protein level) (Combo vs. others)")
p <- p + theme(axis.title.y = element_blank()) + theme(legend.position = "bottom")
p <- p + theme(axis.text.x = element_rect(size = 9))
p

# save output -------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_drug_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
file2write <- paste0(dir_out, "allmodels", "_", treatment.month_tmp, "month", ".pdf")
pdf(file2write, width = 5, height = 2.5, useDingbats = F)
print(p)
dev.off()


