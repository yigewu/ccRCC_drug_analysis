# Yige Wu @ WashU 2021 Mar

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
## input the protein data
protein_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/preprocess/normalize_DIA_protein_data/20210111.v1/RCC_PDX.DIA_Protein.Log2.QuantileNormalized.20210111.v1.tsv")
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")

# preprocess --------------------------------------------------------------
sampleinfo_df <- sampleinfo_df %>%
  mutate(Id_Model = str_split_fixed(string = `Sample ID`, pattern = "_", n = 3)[,1])

# average by model --------------------------------------------------------
avg_pro_df <- protein_df[, c("PG.ProteinAccessions", "PG.Genes", "PG.ProteinDescriptions", "PG.ProteinNames")]
for (id_model_tmp in unique(sampleinfo_df$Id_Model)) {
  ids_sample <- sampleinfo_df$`Sample ID`[sampleinfo_df$Id_Model == id_model_tmp]
  avg_pro_df_tmp <- rowMeans(x = protein_df[, ids_sample], na.rm = T)
  avg_pro_df_tmp <- data.frame(avg_pro_df_tmp)
  colnames(avg_pro_df_tmp) <- id_model_tmp
  ## combine
  avg_pro_df <- cbind(avg_pro_df, avg_pro_df_tmp)
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Average_Protein.", "By_Model.", run_id, ".tsv")
write.table(x = avg_pro_df, file = file2write, quote = F, sep = "\t", row.names = F)
