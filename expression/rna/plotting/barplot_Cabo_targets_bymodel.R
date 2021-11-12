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
foldchanges_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/rna/get_foldchange/get_foldchange_treated_vs_baseline/20211108.v1/RNA_foldchange_treated_vs_baseline.20211108.v1.tsv")
## input detailed sample meta data
meta_data_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/0.sample_info/v3.20210116/RCC_PDX_Samples.20210115.v2.xlsx")

# run by loop -------------------------------------------------------------
meta_data_baseline_df <- meta_data_df %>%
  filter(ShortTag == "Baseline") %>%
  filter(DataType == "RNA-Seq") %>%
  arrange(ModelID)
meta_data_treated_df <- meta_data_df %>%
  filter(grepl(pattern = "Treated", x = ShortTag) | ShortTag == "Control") %>%
  filter(DataType == "RNA-Seq") %>%
  arrange(ModelID)

# prepare plotting parameters ---------------------------------------------
colors_treatment <- c("grey50", RColorBrewer::brewer.pal(n = 5, name = "Set1")[c(1,2,4,3)])
names(colors_treatment) <- c("Baseline", "Treated.Cabo", "Treated.Sap", "Treated.Cabo+Sap", "Control")

# plot by model -----------------------------------------------------------
model_id_tmp <- "RESL5"
treatment_tmp <- "Treated.Cabo"
treatment.month_tmp <- 1
for (model_id_tmp in unique(meta_data_baseline_df$ModelID)) {
  ## prepare plot data
  plotdata_wide_df <- foldchanges_df[, c("genesymbol", paste0(model_id_tmp, "_", c(treatment_tmp, "Control"), "_", treatment.month_tmp, "month"))]
  plotdata_wide_filtered_df <- plotdata_wide_df %>%
    filter(genesymbol %in% genes_rtk_cabo)
  plotdata_wide_filtered_df <- plotdata_wide_filtered_df[rowSums(!is.na(plotdata_wide_filtered_df)) == 3,]
  plotdata_df <- plotdata_wide_filtered_df %>%
    melt() %>%
    filter(!is.infinite(value)) %>%
    mutate(treatment = str_split_fixed(string = variable, pattern = "_", n = 3)[,2]) %>%
    mutate(log2FC = log2(value)) %>%
    arrange(log2FC)
  
  plotdata_df$y_plot <- factor(x = plotdata_df$genesymbol, levels = unique(plotdata_df$genesymbol))
  ## plot
  p <- ggplot()
  p <- p + geom_bar(data = plotdata_df, mapping = aes(x = log2FC, y = y_plot, group = treatment, fill = treatment), stat = "identity", position = 'dodge')
  p <- p + scale_fill_manual(values = colors_treatment)
  p <- p + theme_classic()
  file2write <- paste0(dir_out, model_id_tmp, "_", treatment_tmp, "_", treatment.month_tmp, "month", ".pdf")
  pdf(file2write, width = 7, height = 5.5, useDingbats = F)
  print(p)
  dev.off()
}



