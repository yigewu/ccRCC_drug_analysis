# Yige Wu @ WashU 2021 Nov
## need to overlap with protein changes

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
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
  filter(ShortTag %in% c(treatment_tmp, "Control")) %>%
  filter(DataType == "RNA-Seq") %>%
  arrange(ModelID)

# prepare plotting parameters ---------------------------------------------
colors_treatment <- c("grey50", RColorBrewer::brewer.pal(n = 5, name = "Set1")[c(1,2,4,3)])
names(colors_treatment) <- c("Baseline", "Treated.Cabo", "Treated.Sap", "Treated.Cabo+Sap", "Control")

# plot all model -----------------------------------------------------------
model_id_tmp <- "RESL5"
treatment_tmp <- "Treated.Cabo"
treatment.month_tmp <- 1
model_ids_tmp <- unique(meta_data_baseline_df$ModelID)
model_ids_tmp <- intersect(model_ids_tmp, unique(meta_data_treated_df$ModelID))

## prepare plot data
plotdata_wide_df <- foldchanges_df[, c("genesymbol", 
                                       paste0(model_ids_tmp, "_", treatment_tmp, "_", treatment.month_tmp, "month"),
                                       paste0(model_ids_tmp, "_", "Control", "_", treatment.month_tmp, "month"))]
## filter for sensitivity-related genes
plotdata_wide_filtered_df <- plotdata_wide_df %>%
  filter(RESL5_Treated.Cabo_1month > 0 & !is.infinite(RESL5_Treated.Cabo_1month) & RESL5_Treated.Cabo_1month <= 0.5 & (RESL5_Treated.Cabo_1month/RESL5_Control_1month) <= 0.5)

## filter for resistance-related genes
plotdata_wide_filtered_df <- plotdata_wide_df %>%
  filter(!is.infinite(RESL4_Treated.Cabo_1month) & RESL4_Treated.Cabo_1month >= 2 & (RESL5_Treated.Cabo_1month/RESL5_Control_1month) >= 2)



plotdata_wide_filtered_df <- plotdata_wide_filtered_df[rowSums(!is.na(plotdata_wide_filtered_df)) == ncol(plotdata_wide_filtered_df),]
plotdata_df <- plotdata_wide_filtered_df %>%
  melt() %>%
  filter(!is.infinite(value)) %>%
  mutate(treatment = str_split_fixed(string = variable, pattern = "_", n = 3)[,2]) %>%
  mutate(model_id = str_split_fixed(string = variable, pattern = "_", n = 3)[,1]) %>%
  mutate(log2FC = log2(value))
