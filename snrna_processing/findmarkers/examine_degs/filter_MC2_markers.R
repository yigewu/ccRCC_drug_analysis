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
senescence_df = readxl::read_xlsx(path = "~/Downloads/41467_2022_32552_MOESM4_ESM.xlsx")

# count DEG across samples ------------------------------------------------
deg_count_df = deg_df %>%
  filter(p_val_adj < 0.05) %>%
  select(gene_symbol) %>%
  table() %>%
  as.data.frame() %>%
  rename("gene" = '.') %>%
  arrange(desc(Freq))

write.table(x = deg_count_df, file = paste0(dir_out, "MC2_DEG_count.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)

deg_count_df %>%
  filter(gene %in% senescence_df$`Gene(human)`)
