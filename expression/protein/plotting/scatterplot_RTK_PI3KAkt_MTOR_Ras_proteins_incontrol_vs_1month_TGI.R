# Yige Wu @ WashU 2022 Jan

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
packages = c(
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2",
  "ggrepel",
  "ggpubr"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
source("./ccRCC_drug_analysis/functions.R")
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
tumorvolume_df <- readxl::read_excel(path = "./Resources/Treatment_Lists/Treatment_Summary/Cabo_Sapa_TumorReduction.031522.xlsx", 
                                     sheet = "1month_2month.TGI.v1", 
                                     skip = 2, 
                                     col_names = c("Model",  "Treatment_length", paste0("RTV.", c("vehicle", "Cabo", "Sap", "Cabo+ Sap")), paste0("TGI.", c("Cabo", "Sap", "Cabo+ Sap"))))
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/4.protein/WUSTL to JHU_ccRCC PDX_sample information_01062021_YW.xlsx")

# plot by pathway ---------------------------------------------------------
sampleinfo_df <- sampleinfo_df %>%
  mutate(model_id = str_split_fixed(string = `Sample ID`, pattern = "_", n = 3)[,1]) %>%
  mutate(sample_id = paste0(model_id, "_", Treatment, "_", gsub(x = Treatment_length, pattern = " ", replacement = ""))) %>%
  mutate(group_id = paste0(model_id, "_", gsub(x = Treatment_length, pattern = " ", replacement = ""))) %>%
  filter(Treatment_length == "1 month")
sampleinfo_control_df <- sampleinfo_df %>%
  filter(Treatment == "Con") %>%
  arrange(Treatment, Treatment_length, model_id)
sampleids_control_tmp <- sampleinfo_control_df$`Sample ID`
fontsize_plot <- 16
treatment_tmp <- "Cabo"
for (treatment_tmp in c("Cabo")) {
# for (treatment_tmp in c("Cabo", "Sap", "Cabo+ Sap")) {
  dir_tmp_out <- paste0(dir_out, treatment_tmp, "/")
  dir.create(dir_tmp_out)
  sampleinfo_treated_df <- sampleinfo_df %>%
    filter(Treatment == treatment_tmp) %>%
    arrange(Treatment, Treatment_length, model_id)
  sampleids_treated_tmp <- mapvalues(x = sampleinfo_control_df$model_id, from = sampleinfo_treated_df$model_id, to = as.vector(sampleinfo_treated_df$sample_id))
  
  tumorvolumes_tmp <- mapvalues(x = sampleinfo_control_df$model_id, from = tumorvolume_df$Model, to = unlist(tumorvolume_df[,paste0("TGI.", treatment_tmp)]))
  tumorvolumes_tmp <- as.numeric(tumorvolumes_tmp)
  # for (pathwayname_tmp in c("RTK_ligand")) {
    for (pathwayname_tmp in c("RTK_ligand", "PI3K_Akt", "Ras_MAPK")) {
    pathway2members_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Gene_Lists/Pathway_score_members.011222.xlsx", sheet = pathwayname_tmp)
    pathway2members_df <- pathway2members_df %>%
      mutate(ID = paste0(Gene_symbol, "_", ifelse(Is_phospho == "Yes", Site_phospho, "Protein"))) %>%
      filter(ID %in% exp_df$ID)
    
    # for (protein_tmp in "MET_Protein") {
    for (protein_tmp in pathway2members_df$ID) {
      exp_tmp <- unlist(exp_df[exp_df$ID == protein_tmp, sampleids_control_tmp])
      
      ## make plot data
      plotdata_df <- data.frame(sample_id = sampleids_control_tmp,
                                model_id = str_split_fixed(string = sampleids_treated_tmp, pattern = "_", n = 3)[,1],
                                treatment_length = str_split_fixed(string = sampleids_treated_tmp, pattern = "_", n = 3)[,3],
                                exp_value = exp_tmp, 
                                relative_tumor_volume = tumorvolumes_tmp)
      plotdata_df <- plotdata_df %>%
        mutate(x_plot = exp_value) %>%
        mutate(y_plot = relative_tumor_volume) %>%
        filter(!is.na(x_plot) & !is.na(y_plot))

      plotdata_tmp_df <- plotdata_df %>%
        filter(treatment_length == "1month") 
      p <- ggscatter(data = plotdata_tmp_df, x = "exp_value", y = "relative_tumor_volume",
                     add = "reg.line",  # Add regressin line
                     label = "model_id", font.label = c(fontsize_plot, "plain"),
                     add.params = list(color = "grey", fill = "lightgray", linetype = 2), # Customize reg. line
                     conf.int = T # Add confidence interval
      )
      p <- p + stat_cor(method = "pearson",
                        label.x = min(plotdata_tmp_df$x_plot),
                        label.y = (max(plotdata_tmp_df$y_plot) + 0.1), size = 7)
      # p <- p + ggtitle(label = paste0(protein_tmp, "~", treatment_tmp, " 1-month response"))
      p <- p + theme_classic(base_size = fontsize_plot)
      p <- p + xlab(paste0(gsub(x = protein_tmp, pattern = "_", replacement = " "), " level (baseline)"))
      p <- p + ylab(paste0("Tumor growth inhibition (day 28)"))
      p <- p + xlim(c(min(plotdata_tmp_df$x_plot)-0.11, max(plotdata_tmp_df$x_plot)+0.11))
      p <- p + ylim(c(min(plotdata_tmp_df$y_plot)-0.05, max(plotdata_tmp_df$y_plot)+0.1))
      p <- p + theme(axis.text = element_text(color = "black", size = fontsize_plot))
      # file2write <- paste0(dir_tmp_out, protein_tmp, ".", pathwayname_tmp, ".", treatment_tmp, ".", "1month", ".png")
      # png(file2write, width = 600, height = 500, res = 150)
      # print(p)
      # dev.off()
      file2write <- paste0(dir_tmp_out, protein_tmp, ".", pathwayname_tmp, ".", treatment_tmp, ".", "1month", ".pdf")
      pdf(file2write, width = 4, height = 4, useDingbats = F)
      print(p)
      dev.off()
      
    }
  }
}
