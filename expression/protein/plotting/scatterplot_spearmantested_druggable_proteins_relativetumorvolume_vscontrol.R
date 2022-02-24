# Yige Wu @ WashU 2022 Jan

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
source("./ccRCC_drug_analysis/plotting.R")
library(ggpubr)
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
enricher_top_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/pathway/unite_ora_msigdb_H_CP_treated_relativetumorvolume_human_protein_markers_incontrol/20220203.v1/ora_msigdb_H_CP_human_protein_markers_in_control.top.20220203.v1.tsv")
druggable_df <- fread(data.table = F, input = "./Resources/Knowledge/Gene_Lists/approved_target_ids_all.csv")

# plot by pathway ---------------------------------------------------------
sampleinfo_df <- sampleinfo_df %>%
  mutate(model_id = str_split_fixed(string = `Sample ID`, pattern = "_", n = 3)[,1]) %>%
  mutate(sample_id = paste0(model_id, "_", Treatment, "_", gsub(x = Treatment_length, pattern = " ", replacement = ""))) %>%
  mutate(group_id = paste0(model_id, "_", gsub(x = Treatment_length, pattern = " ", replacement = "")))
sampleinfo_control_df <- sampleinfo_df %>%
  filter(Treatment == "Con") %>%
  arrange(Treatment, Treatment_length, model_id)
sampleids_control_tmp <- sampleinfo_control_df$`Sample ID`
enricher_top_df <- enricher_top_df %>%
  filter(Count >=3) %>%
  mutate(treatment = str_split_fixed(string = test, pattern = "_", n = 9)[,5]) %>%
  mutate(treatment_matched = ifelse(treatment == "CaboSap", "Cabo+ Sap", treatment)) %>%
  mutate(assoc_direction = str_split_fixed(string = test, pattern = "_", n = 9)[,8])

pathwayname_tmp <- "Cell_cycle"
treatment_tmp <- "Sap"
# for (treatment_tmp in c("Sap")) {
for (treatment_tmp in c("Cabo", "Cabo+ Sap", "Sap")) {
  dir_tmp_out <- paste0(dir_out, treatment_tmp, "/")
  dir.create(dir_tmp_out)
  sampleinfo_treated_df <- sampleinfo_df %>%
    filter(Treatment == treatment_tmp) %>%
    arrange(Treatment, Treatment_length, model_id)
  sampleids_treated_tmp <- mapvalues(x = sampleinfo_control_df$group_id, from = sampleinfo_treated_df$group_id, to = as.vector(sampleinfo_treated_df$sample_id))
  
  tumorvolumes_tmp <- mapvalues(x = sampleids_treated_tmp, from = paste0(tumorvolume_df$Model, "_", treatment_tmp, "_", tumorvolume_df$Treatment_length), to = unlist(tumorvolume_df[,paste0(treatment_tmp, " vs. Control")]))
  tumorvolumes_tmp <- as.numeric(tumorvolumes_tmp)
  
  ## get genes to highlight
  genestrings_tmp <- enricher_top_df$geneID[enricher_top_df$treatment_matched == treatment_tmp]
  genes_list_tmp <- str_split(string = genestrings_tmp, pattern = "\\/")
  genes_tmp <- unlist(genes_list_tmp); genes_tmp <- unique(genes_tmp)
  genes_tmp <- genes_tmp[genes_tmp %in% druggable_df$`Gene Name`]
  
  for (protein_tmp in paste0(genes_tmp, "_Protein")) {
    # for (protein_tmp in pathway2members_df$ID) {
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
    # ## plot
    # p <- ggscatter(data = plotdata_df, x = "exp_value", y = "relative_tumor_volume", shape = "treatment_length",
    #                add = "reg.line",  # Add regressin line
    #                # add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    #                conf.int = TRUE # Add confidence interval
    # )
    # p <- p + stat_cor(method = "pearson",
    #                   label.x = quantile(x = plotdata_df$exp_value, probs = 0.001, na.rm = T))
    # p <- p + ggtitle(label = paste0(protein_tmp, "~", treatment_tmp, " response"))
    # p <- p + xlab(paste0(protein_tmp, " abundance in control"))
    # p <- p + ylab(paste0("Relative tumor volume (", treatment_tmp, "/Control)"))
    # file2write <- paste0(dir_tmp_out, protein_tmp, ".", pathwayname_tmp, ".", treatment_tmp, ".png")
    # png(file2write, width = 600, height = 600, res = 150)
    # print(p)
    # dev.off()
    # 
    # p <- ggscatter(data = plotdata_df, x = "exp_value", y = "relative_tumor_volume", shape = "treatment_length", color = "treatment_length",
    #                add = "reg.line",  # Add regressin line
    #                # add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    #                conf.int = TRUE # Add confidence interval
    # )
    # p <- p + stat_cor(aes(color = treatment_length),
    #                   method = "pearson",
    #                   label.x = quantile(x = plotdata_df$exp_value, probs = 0.001, na.rm = T))
    # p <- p + ggtitle(label = paste0(protein_tmp, "~", treatment_tmp, " response"))
    # p <- p + xlab(paste0(protein_tmp, " abundance in control"))
    # p <- p + ylab(paste0("Relative tumor volume (", treatment_tmp, "/Control)"))
    # file2write <- paste0(dir_tmp_out, protein_tmp, ".", treatment_tmp, ".bytreatmentlength.png")
    # png(file2write, width = 600, height = 600, res = 150)
    # print(p)
    # dev.off()
    
    plotdata_tmp_df <- plotdata_df %>%
      filter(treatment_length == "1month") 
    p <- ggscatter(data = plotdata_tmp_df, x = "exp_value", y = "relative_tumor_volume",
                   add = "reg.line",  # Add regressin line
                   label = "model_id",
                   add.params = list(color = "grey", fill = "lightgray", linetype = 2), # Customize reg. line
                   conf.int = T # Add confidence interval
    )
    p <- p + stat_cor(method = "pearson",
                      label.x = min(plotdata_tmp_df$x_plot),
                      label.y = (max(plotdata_tmp_df$y_plot) + 0.1))
    p <- p + ggtitle(label = paste0(protein_tmp, "~", treatment_tmp, " 1-month response"))
    p <- p + xlab(paste0(protein_tmp, " abundance in control"))
    p <- p + ylab(paste0("Relative tumor volume (", treatment_tmp, "/Control)"))
    p <- p + xlim(c(min(plotdata_tmp_df$x_plot)-0.11, max(plotdata_tmp_df$x_plot)+0.11))
    p <- p + ylim(c(min(plotdata_tmp_df$y_plot)-0.05, max(plotdata_tmp_df$y_plot)+0.1))
    
    file2write <- paste0(dir_tmp_out, protein_tmp, ".", treatment_tmp, ".", "1month", ".png")
    png(file2write, width = 600, height = 500, res = 150)
    print(p)
    dev.off()
    
  }
  
}
