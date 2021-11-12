# Yige Wu @ WashU 2021 Nov

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
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
## input fold changes
foldchanges_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/rna/get_foldchange/get_foldchange_treated_vs_baseline/20211109.v1/RNA_foldchange_treated_vs_baseline.20211109.v1.tsv")
## input detailed sample meta data
meta_data_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/0.sample_info/v3.20210116/RCC_PDX_Samples.20210115.v2.xlsx")
## input the markers
# markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/rna/get_foldchange/filter_Sap_sensitive_resistant_genes_overlapping_proteins/20211110.v1/Sap_related_genes.20211110.v1.tsv")
markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/rna/get_foldchange/filter_Sap_sensitive_resistant_genes_overlapping_proteins/20211110.v3/Sap_related_genes.20211110.v3.tsv")

# preprocess -------------------------------------------------------------
genes_process <- markers_df$genesymbol[markers_df$gene_category == "Sap resistant"]
meta_data_baseline_df <- meta_data_df %>%
  filter(ShortTag == "Baseline") %>%
  filter(DataType == "RNA-Seq") %>%
  arrange(ModelID)
meta_data_treated_df <- meta_data_df %>%
  filter(ShortTag %in% c(treatment_tmp, "Control")) %>%
  filter(DataType == "RNA-Seq") %>%
  arrange(ModelID)

# prepare plotting parameters ---------------------------------------------
colors_treatment <- c("grey50", RColorBrewer::brewer.pal(n = 5, name = "Set1")[c(1,2,4,3)])
names(colors_treatment) <- c("Baseline", "Treated.Cabo", "Treated.Sap", "Treated.Cabo+Sap", "Control")
cap_value <- 4

# plot all model -----------------------------------------------------------
treatment_tmp <- "Treated.Sap"
treatment.month_tmp <- 1
model_ids_tmp <- unique(meta_data_baseline_df$ModelID)
model_ids_tmp <- intersect(model_ids_tmp, unique(meta_data_treated_df$ModelID))

## prepare plot data
# plotdata_wide_df <- foldchanges_df[, c("genesymbol", 
#                                        paste0(model_ids_tmp, "_", treatment_tmp, "_", treatment.month_tmp, "month"),
#                                        paste0(model_ids_tmp, "_", "Control", "_", treatment.month_tmp, "month"))]
plotdata_wide_df <- foldchanges_df[, c("genesymbol", 
                                       paste0(model_ids_tmp, "_", "Treated.Sap", "_", treatment.month_tmp, "month"),
                                       paste0(model_ids_tmp, "_", "Treated.Cabo", "_", treatment.month_tmp, "month"),
                                       paste0(model_ids_tmp, "_", "Control", "_", treatment.month_tmp, "month"))]

plotdata_wide_filtered_df <- plotdata_wide_df %>%
  filter(genesymbol %in% genes_process)
plotdata_wide_filtered_df <- plotdata_wide_filtered_df[rowSums(!is.na(plotdata_wide_filtered_df)) == ncol(plotdata_wide_filtered_df),]
plotdata_df <- plotdata_wide_filtered_df %>%
  melt() %>%
  filter(!is.infinite(value)) %>%
  mutate(treatment = str_split_fixed(string = variable, pattern = "_", n = 3)[,2]) %>%
  mutate(model_id = str_split_fixed(string = variable, pattern = "_", n = 3)[,1]) %>%
  mutate(log2FC = log2(value)) %>%
  mutate(x_plot = ifelse(log2FC > cap_value, cap_value,
                         ifelse(log2FC < -cap_value, -cap_value, log2FC)))
orderdata_df <- plotdata_df %>%
  filter(model_id == "RESL10") %>%
  filter(treatment == "Treated.Sap") %>%
  arrange(log2FC)
plotdata_df$model_id <- factor(x = plotdata_df$model_id, levels = c("RESL5", "RESL11", "RESL3", "RESL4", "RESL10"))
plotdata_df$y_plot <- factor(x = plotdata_df$genesymbol, levels = unique(orderdata_df$genesymbol))

## plot
p <- ggplot(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot))
p <- p + geom_bar(mapping = aes(fill = treatment), stat = "identity", position = 'dodge')
p <- p + scale_fill_manual(values = colors_treatment)
p <- p + facet_grid(cols = vars(model_id))
p <- p + theme_classic()
p <- p + theme(axis.title.y = element_blank()) + theme(legend.position = "bottom")
p <- p + xlim(c(-cap_value, cap_value))
p <- p + xlab("Log2(Expression-treated vs. Expression-baseline)")
p <- p + theme(axis.title.y = element_blank()) + theme(legend.position = "bottom")

file2write <- paste0(dir_out, "allmodels", "_", treatment_tmp, "_", treatment.month_tmp, "month", ".pdf")
pdf(file2write, width = 5, height = 3.5, useDingbats = F)
print(p)
dev.off()



