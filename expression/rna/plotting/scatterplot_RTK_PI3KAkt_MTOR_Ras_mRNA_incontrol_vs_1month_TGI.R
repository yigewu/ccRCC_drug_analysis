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
exp_df <- fread("./Data_Freeze/v1.dataFreeze.washU_rcc/3.geneExp/v3.20210116/datafreeze.v3.kallisto.geneExp.protein_coding.tsv", data.table = F)
## input the relative tumor volume
tumorvolume_df <- readxl::read_excel(path = "./Resources/Treatment_Lists/Treatment_Summary/Cabo_Sapa_TumorReduction.031522.xlsx", 
                                     sheet = "1month_2month.TGI.v1", 
                                     skip = 2, 
                                     col_names = c("Model",  "Treatment_length", paste0("RTV.", c("vehicle", "Cabo", "Sap", "Cabo+Sap")), paste0("TGI.", c("Cabo", "Sap", "Cabo+Sap"))))
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/0.sample_info/v3.20210116/RCC_PDX_Samples.20210115.v2.xlsx")

# plot by pathway ---------------------------------------------------------
sampleinfo_df <- sampleinfo_df %>%
  mutate(Treatment = ShortTag)
sampleinfo_control_df <- sampleinfo_df %>%
  filter(Treatment == "Control") %>%
  filter(Treatment.Month == 1) %>%
  filter(DataType == "RNA-Seq")
analysisids_control_tmp <- sampleinfo_control_df$Analysis_ID
tumorvolume_df <- tumorvolume_df %>%
  filter(Treatment_length == 28)
fontsize_plot <- 16

treatment_tmp <- "Cabo"
for (treatment_tmp in c("Cabo")) {
# for (treatment_tmp in c("Cabo", "Sap", "Cabo+Sap")) {
  dir_tmp_out <- paste0(dir_out, treatment_tmp, "/")
  dir.create(dir_tmp_out)
  tumorvolumes_tmp <- as.numeric(unlist(tumorvolume_df[, paste0("TGI.", treatment_tmp)]))
  tumorvolumes_df <- data.frame(model_id = tumorvolume_df$Model,
                                relative_tumor_volume = tumorvolumes_tmp)
  
  # for (pathwayname_tmp in c("RTK_ligand")) {
    for (pathwayname_tmp in c("RTK_ligand", "PI3K_Akt", "Ras_MAPK")) {
    pathway2members_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Gene_Lists/Pathway_score_members.011222.xlsx", sheet = pathwayname_tmp)
    pathway2members_df <- pathway2members_df %>%
      mutate(ID = paste0(Gene_symbol, "_", ifelse(Is_phospho == "Yes", Site_phospho, "Protein"))) %>%
      filter(Gene_symbol %in% exp_df$Name)
    
    # for (protein_tmp in "MET_Protein") {
    for (gene_tmp in unique(pathway2members_df$Gene_symbol)) {
      exp_tmp <- unlist(exp_df[exp_df$Name == gene_tmp, analysisids_control_tmp])
      exp_tmp_df <- sampleinfo_control_df
      exp_tmp_df$exp_tmp_value <- log2(exp_tmp+1)
      exp_avg_df2 <- exp_tmp_df %>%
        filter(Treatment.Month == 1) %>%
        group_by(ModelID) %>%
        summarise(exp_value = mean(exp_tmp_value, na.rm = T))
      
      ## make plot data
      plotdata_df <- merge(x = tumorvolumes_df, y = exp_avg_df2, by.x = c("model_id"), by.y = c("ModelID"), all.x = T)
      plotdata_df <- plotdata_df %>%
        mutate(x_plot = exp_value) %>%
        mutate(y_plot = relative_tumor_volume) %>%
        filter(!is.na(x_plot) & !is.na(y_plot))

      p <- ggscatter(data = plotdata_df, x = "exp_value", y = "relative_tumor_volume",
                     add = "reg.line",  # Add regressin line
                     label = "model_id", font.label = c(fontsize_plot, "plain"),
                     add.params = list(color = "grey", fill = "lightgray", linetype = 2), # Customize reg. line
                     conf.int = T # Add confidence interval
      )
      p <- p + stat_cor(method = "pearson",
                        label.x = min(plotdata_df$x_plot),
                        label.y = (max(plotdata_df$y_plot) + 0.1), size = 7)
      # p <- p + ggtitle(label = paste0(protein_tmp, "~", treatment_tmp, " 1-month response"))
      p <- p + theme_classic(base_size = fontsize_plot)
      p <- p + xlab(paste0(gene_tmp, " gene expression (baseline)"))
      p <- p + ylab(paste0("Tumor growth inhibition (day 28)"))
      p <- p + xlim(c(min(plotdata_df$x_plot)-0.11, max(plotdata_df$x_plot)+0.11))
      p <- p + ylim(c(min(plotdata_df$y_plot)-0.05, max(plotdata_df$y_plot)+0.1))
      p <- p + theme(axis.text = element_text(color = "black", size = fontsize_plot))
      # file2write <- paste0(dir_tmp_out, protein_tmp, ".", pathwayname_tmp, ".", treatment_tmp, ".", "1month", ".png")
      # png(file2write, width = 600, height = 500, res = 150)
      # print(p)
      # dev.off()
      file2write <- paste0(dir_tmp_out, gene_tmp, ".", pathwayname_tmp, ".", treatment_tmp, ".", "1month", ".pdf")
      pdf(file2write, width = 4, height = 4, useDingbats = F)
      print(p)
      dev.off()
      
    }
  }
}
