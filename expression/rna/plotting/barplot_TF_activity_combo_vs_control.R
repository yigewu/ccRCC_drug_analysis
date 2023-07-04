# Yige Wu @ WashU 2023 Jun

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
packages = c(
  "plyr",
  "stringr",
  "reshape2",
  "data.table",
  "dplyr",
  "ggpubr",
  "rstatix"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

# input data --------------------------------------------------------------
## input TF activty per sample
# exp_df = fread(data.table = F, input = "./Resources/Analysis_Results/expression/rna/pathway/calculate_TF_activity_bysample/20230621.v2/TF_activity.long.20230621.v2.tsv")
exp_df = fread(data.table = F, input = "./Resources/Analysis_Results/expression/rna/pathway/calculate_TF_activity_bysample/20230621.v3/TF_activity.long.20230621.v3.tsv")
## input TF activity test
exp_test_df = fread(data.table = F, input = "./Resources/Analysis_Results/expression/rna/test/ttest_paired_diff_TFactivity_1month_combo_vs_single/20230621.v4/TFactivity.Ttest.Paired.1month.Combo_vs_single.20230621.v4.tsv")

## input ERK downstream
genes2filter_df = readxl::read_xlsx(path = "./Resources/Knowledge/41580_2020_255_MOESM1_ESM.xlsx", skip = 1)

# select TFs to plot ------------------------------------------------------
length(unique(exp_df$TF_genesymbol))
genes_select_df = exp_test_df %>%
  filter(group2 == "Control") %>%
  filter(Name %in% genes2filter_df$SUBSTRATE_GENE) %>%
  filter(diff_estimate < 0) %>% ## Combo < Control
  filter(!(Name %in% c("FOXO3", "SMAD4", "POU5F1", "MITF", "STAT3"))) %>%
  arrange(pvalue)

genes_filter1_df = exp_test_df %>%
  filter(group2 == "Treated.Cabo") %>%
  filter(Name %in% genes2filter_df$SUBSTRATE_GENE) %>%
  filter(diff_estimate < 0)

genes_filter2_df = exp_test_df %>%
  filter(group2 == "Treated.Sap") %>%
  filter(Name %in% genes2filter_df$SUBSTRATE_GENE) %>%
  filter(diff_estimate < 0)

genes_select_df = genes_select_df %>%
  filter(Name %in% genes_filter1_df$Name) %>%
  filter(Name %in% genes_filter2_df$Name) %>%
  filter(pvalue < 0.2)

genes_select_df$Name

# make plot data ----------------------------------------------------------
plot_data_df = exp_df %>%
  filter(TF_genesymbol %in% genes_select_df$Name) %>%
  filter(!(PairTag %in% c("RESL5D.27", "RESL5E.28", "RESL5E.30", "RESL10F.38"))) %>%
  filter(ShortTag %in% c("Control", "Treated.Cabo+Sap", "Treated.Cabo", "Treated.Sap")) %>%
  filter(Treatment.Month == 1) %>%
  mutate(Treatment_length = paste0(Treatment.Month, " month")) %>%
  mutate(PairTag.middle = str_split_fixed(string = analysis_id, pattern = "_", n = 3)[,2]) %>%
  mutate(PairTag = ifelse(!is.na(PairTag), PairTag, paste0(ModelID, "_", PairTag.middle))) %>%
  mutate(id_model_length = paste0(PairTag, "_", Treatment_length))
plot_data_df$Treatment = mapvalues(plot_data_df$ShortTag, 
                                   from = c("Control", "Treated.Cabo", "Treated.Sap", "Treated.Cabo+Sap"),
                                   to = c("Control", "Cabozantinib", "Sapanisertib", "Cabozantinib+Sapanisertib"))
plot_data_df$Treatment = factor(x = plot_data_df$Treatment, levels = c("Control", "Cabozantinib", "Sapanisertib", "Cabozantinib+Sapanisertib"))
ids_proces = unique(plot_data_df$id_model_length[plot_data_df$Treatment == "Cabozantinib+Sapanisertib"])
plot_data_df = plot_data_df  %>%
  filter(id_model_length %in% ids_proces) %>%
  arrange(TF_genesymbol, id_model_length, desc(Treatment))
unique(plot_data_df$Treatment)
plot_data_df$TF_genesymbol = factor(plot_data_df$TF_genesymbol, levels = genes_select_df$Name)
# plot_data_df$y_plot <- (plot_data_df$TF_activity_score - min(plot_data_df$TF_activity_score))/(max(plot_data_df$TF_activity_score) - min(plot_data_df$TF_activity_score))
# plot_data_df$y_plot <- (plot_data_df$TF_activity_score - min(plot_data_df$TF_activity_score))
plot_data_df$y_plot <- plot_data_df$TF_activity_score

stat.test <- plot_data_df %>%
  group_by(TF_genesymbol) %>%
  t_test(y_plot ~ Treatment, paired = T, 
         comparisons = list(c("Cabozantinib+Sapanisertib", "Control"))) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p")
stat.test <- stat.test %>% add_xy_position(fun = "mean", x = "Treatment")
stat.test$y.position[stat.test$TF_genesymbol == "ETV1"] = 0.05

## make colors
color_red <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[1]
color_green <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[3]
color_yellow <- RColorBrewer::brewer.pal(n = 7, name = "Set2")[6]
color_grey <- "grey50"
colors_plot <- c("Control" = color_grey, "Cabozantinib" = color_red, "Sapanisertib" = color_green, "Cabozantinib+Sapanisertib" = color_yellow)

# plot --------------------------------------------------------------------
p <- ggbarplot(data = plot_data_df, 
               x = "Treatment", y = "y_plot", fill = "Treatment",
               # add = c("mean", "dotplot"),
               add = c("mean"),
               position = position_dodge(),
               facet.by = "TF_genesymbol", nrow = 2)
p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01)
p <- p + geom_hline(yintercept = 0)
p <- p + ylim(c(-0.3, 0.2))
p <- p + scale_fill_manual(values = colors_plot)
p <- p + labs(y = "Transcription factor activity score")
p <- p + guides(fill = guide_legend(nrow = 2))
p <- p + theme_classic(base_size = 15)
p <- p + theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.line.x = element_blank(),
               axis.title.x = element_blank(),
               legend.position = "bottom")
# save output -------------------------------------------------------------
## set run id
version_tmp <- 3
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_drug_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
file2write = paste0(dir_out, "14TFs.pdf")
pdf(file2write, width = 6, height = 5.5, useDingbats = F)
print(p)
dev.off()
