# Yige Wu @WashU Apr 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input human expression
exp_hs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/average_expression/avgexp_sct_data_humancells_on_katmai/20210406.v1/HumanCells.AverageExpression.20210406.v1.tsv")
## input mouse expression
exp_mm_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/average_expression/avgexp_sct_data_mousecells_bycelltypeshort_on_katmai/20210406.v1/AverageExpression.ByCellTypeShorter.20210406.v1.tsv")

# merge by human ortholog -------------------------------------------------
exp_hs_df <- exp_hs_df %>%
  rename(genesymbol_human = V1)
exp_mm_df <- exp_mm_df %>%
  rename(genesymbol_mouse = V1) %>%
  mutate(genesymbol_human = toupper(genesymbol_mouse))
colnames(exp_mm_df) <- gsub(x = colnames(exp_mm_df), pattern = "SCT\\.", replacement = "")
exp_merged_df <- merge(x = exp_hs_df, 
                       y = exp_mm_df,
                       by = c("genesymbol_human"), all = T)
exp_merged_df <- exp_merged_df %>%
  select(genesymbol_human, genesymbol_mouse, Tumor_cells, Endothelial.cells, Myofibroblasts, Macrophages, Fibroblasts)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "SCT.data.", "AverageExpression.", run_id, ".tsv")
write.table(x = exp_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)
