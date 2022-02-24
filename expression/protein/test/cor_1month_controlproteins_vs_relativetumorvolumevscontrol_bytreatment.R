# Yige Wu @ WashU 2022 Jan

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
library(doParallel)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input pathway score changes for treated vs. control
exp_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/unite_protein_phospho_bylog2intensity/20220120.v1/Protein_Phospho_Log2Intensity.20220120.v1.tsv")
## input the relative tumor volume
tumorvolume_df <- readxl::read_excel(path = "./Resources/Treatment_Lists/Treatment_Summary/Cabo_Sapa_TumorReduction.012022.xlsx", sheet = "1month_2month_divide")
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")

# preprocess ---------------------------------------------------------
## reformat sample info
sampleinfo_df <- sampleinfo_df %>%
  mutate(model_id = str_split_fixed(string = `Sample ID`, pattern = "_", n = 3)[,1]) %>%
  mutate(sample_id = paste0(model_id, "_", Treatment, "_", gsub(x = Treatment_length, pattern = " ", replacement = ""))) %>%
  mutate(group_id = paste0(model_id, "_", gsub(x = Treatment_length, pattern = " ", replacement = "")))
sampleinfo_control_df <- sampleinfo_df %>%
  filter(Treatment == "Con") %>%
  filter(Treatment_length == "1 month") %>%
  arrange(Treatment, Treatment_length, model_id)
sampleids_control_tmp <- sampleinfo_control_df$`Sample ID`
## reformat expression
nrow(exp_df)
datapoints_cutoff <- 5
exp_process_df <- exp_df[rowSums(!is.na(exp_df[,sampleids_control_tmp])) >= datapoints_cutoff,]
nrow(exp_process_df)

# test --------------------------------------------------------------------
treatment_tmp <- "Cabo"
no_cores <- 4
registerDoParallel(cores = no_cores)
cor_result_sup_df <- NULL
# for (treatment_tmp in c("Cabo")) {
for (treatment_tmp in c("Cabo", "Sap", "Cabo+ Sap")) {
  sampleinfo_treated_df <- sampleinfo_df %>%
    filter(Treatment == treatment_tmp) %>%
    filter(Treatment_length == "1 month") %>% 
    arrange(Treatment, Treatment_length, model_id)
  sampleids_treated_tmp <- mapvalues(x = sampleinfo_control_df$group_id, from = sampleinfo_treated_df$group_id, to = as.vector(sampleinfo_treated_df$sample_id))
  
  tumorvolumes_tmp <- mapvalues(x = sampleids_treated_tmp, from = paste0(tumorvolume_df$Model, "_", treatment_tmp, "_", tumorvolume_df$Treatment_length), to = unlist(tumorvolume_df[,paste0(treatment_tmp, " vs. Control")]))
  tumorvolumes_tmp <- as.numeric(tumorvolumes_tmp)
  
  # start_time <- Sys.time()
  cor_result_list<- foreach(protein_tmp = exp_process_df$ID) %dopar% {
    # for (protein_tmp in pathway2members_df$ID) {
    exp_tmp <- unlist(exp_process_df[exp_process_df$ID == protein_tmp, sampleids_control_tmp])
    testdata_df <- data.frame(exp_value = exp_tmp, 
                              relative_tumor_volume = tumorvolumes_tmp)
    testdata_df <- testdata_df %>%
      filter(!is.na(exp_value))
    result_obj <- cor.test(x = testdata_df$exp_value, y = testdata_df$relative_tumor_volume, method = "spearman")
    result_vec <- c(result_obj$p.value, result_obj$estimate)
    return(result_vec)
  }
  # end_time <- Sys.time()
  # end_time - start_time
  cor_result_df <- do.call(rbind.data.frame, cor_result_list)
  colnames(cor_result_df) <- c("p.value", "rho")
  cor_result_df$ID <- exp_process_df$ID
  cor_result_df$treatment = treatment_tmp
  cor_result_sup_df <- rbind(cor_result_sup_df, cor_result_df)
}

# post-test processing ----------------------------------------------------
cor_result_sup_df <- cor_result_sup_df %>%
  mutate(abundance_type = ifelse(grepl(pattern = "Protein", x = ID), "total_protein", "PTM"))
cor_result_sup_df <- unique(cor_result_sup_df)
cor_result_sup_df$fdr_byabundancetype[cor_result_sup_df$abundance_type == "total_protein"] <- p.adjust(cor_result_sup_df$p.value[cor_result_sup_df$abundance_type == "total_protein"])
cor_result_sup_df$fdr_byabundancetype[cor_result_sup_df$abundance_type != "total_protein"] <- p.adjust(cor_result_sup_df$p.value[cor_result_sup_df$abundance_type != "total_protein"])

# write file --------------------------------------------------------------
file2write <- paste0(dir_out, "spearman.1month.controlproteins_vs_relativetumorvolume.", run_id, ".tsv")
write.table(x = cor_result_sup_df, file = file2write, quote = F, sep = "\t", row.names = F)
