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
## input expression
exp_df <- fread("./Data_Freeze/v1.dataFreeze.washU_rcc/3.geneExp/v3.20210116/datafreeze.v3.kallisto.geneExp.protein_coding.tsv", data.table = F)
## input the relative tumor volume
tumorvolume_df <- readxl::read_excel(path = "./Resources/Treatment_Lists/Treatment_Summary/Cabo_Sapa_TumorReduction.012022.xlsx", sheet = "1month_2month_divide")
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/0.sample_info/v3.20210116/RCC_PDX_Samples.20210115.v2.xlsx")

# plot by pathway ---------------------------------------------------------
baseline_tmp <- "Control"
treatment_tmp <- "Cabo"
gene_tmp <- "MET"
# sampleinfo_df <- sampleinfo_df %>%
#   mutate(sample_id = paste0(model_id, "_", Treatment, "_", gsub(x = Treatment_length, pattern = " ", replacement = ""))) %>%
#   mutate(group_id = paste0(model_id, "_", gsub(x = Treatment_length, pattern = " ", replacement = "")))
sampleinfo_baseline_df <- sampleinfo_df %>%
  filter(ShortTag == baseline_tmp) %>%
  filter(DataType == "RNA-Seq")
analysisids_baseline_tmp <- sampleinfo_baseline_df$Analysis_ID
# for (treatment_tmp in c("Sap")) {
for (treatment_tmp in c("Cabo", "Sap", "Cabo+ Sap")) {
  dir_tmp_out <- paste0(dir_out, treatment_tmp, "/")
  dir.create(dir_tmp_out)
  tumorvolumes_tmp <- as.numeric(unlist(tumorvolume_df[,paste0(treatment_tmp, " vs. Control")]))
  tumorvolumes_df <- data.frame(model_id = tumorvolume_df$Model,
                            treatment_length = tumorvolume_df$Treatment_length,
                            relative_tumor_volume = tumorvolumes_tmp)
  
  # for (pathwayname_tmp in c("TSC_mTOR")) {
    for (pathwayname_tmp in c("RTK_ligand", "TSC_mTOR", "PI3K_Akt", "Ras_MAPK")) {
    pathway2members_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Gene_Lists/Pathway_score_members.011222.xlsx", sheet = pathwayname_tmp)
    pathway2members_df <- pathway2members_df %>%
      mutate(ID = paste0(Gene_symbol, "_", ifelse(Is_phospho == "Yes", Site_phospho, "Protein"))) %>%
      # filter(Is_phospho == "No") %>%
      filter(Gene_symbol %in% exp_df$Name)
    
    # for (gene_tmp in "GLUD2") {
    for (gene_tmp in unique(pathway2members_df$Gene_symbol)) {
      exp_tmp <- unlist(exp_df[exp_df$Name == gene_tmp, analysisids_baseline_tmp])
      exp_tmp_df <- sampleinfo_baseline_df
      exp_tmp_df$exp_tmp_value <- log2(exp_tmp+1)
      exp_avg_df2 <- exp_tmp_df %>%
        mutate(Treatment_length = paste0(Treatment.Month, "month")) %>%
        group_by(ModelID, Treatment_length) %>%
        summarise(exp_value = mean(exp_tmp_value, na.rm = T))
      ## make plot data
      plotdata_df <- merge(x = tumorvolumes_df, y = exp_avg_df2, by.x = c("model_id", "treatment_length"), by.y = c("ModelID", "Treatment_length"), all.x = T)
      plotdata_df <- plotdata_df %>%
        mutate(x_plot = exp_value) %>%
        mutate(y_plot = relative_tumor_volume) %>%
        filter(!is.na(x_plot) & !is.na(y_plot))
      ## plot
      p <- ggscatter(data = plotdata_df, x = "exp_value", y = "relative_tumor_volume", shape = "treatment_length",
                     add = "reg.line",  # Add regressin line
                     # add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                     conf.int = F # Add confidence interval
      )
      p <- p + stat_cor(method = "pearson",
                        label.x = min(plotdata_tmp_df$x_plot, na.rm = T),
                        label.y = (max(plotdata_tmp_df$y_plot, na.rm = T) + 0.02))
      p <- p + ggtitle(label = paste0(gene_tmp, " mRNA ~", treatment_tmp, " response"))
      p <- p + xlab(paste0(gene_tmp, " gene expression at baseline"))
      p <- p + ylab(paste0("Relative tumor volume (", treatment_tmp, "/Control)"))
      file2write <- paste0(dir_tmp_out, gene_tmp, ".", pathwayname_tmp, ".", treatment_tmp, ".png")
      png(file2write, width = 600, height = 600, res = 150)
      print(p)
      dev.off()
      
      p <- ggscatter(data = plotdata_df, x = "exp_value", y = "relative_tumor_volume", shape = "treatment_length", color = "treatment_length",
                     add = "reg.line",  # Add regressin line
                     add.params = list(linetype = 2), # Customize reg. line
                     conf.int = F # Add confidence interval
      )
      p <- p + stat_cor(aes(color = treatment_length),
                        method = "pearson")
      p <- p + ggtitle(label = paste0(gene_tmp, " mRNA ~", treatment_tmp, " response"))
      p <- p + xlab(paste0(gene_tmp, " gene expression at baseline"))
      p <- p + ylab(paste0("Relative tumor volume (", treatment_tmp, "/Control)"))
      file2write <- paste0(dir_tmp_out, gene_tmp, ".", pathwayname_tmp, ".", treatment_tmp, ".bytreatmentlength.png")
      png(file2write, width = 600, height = 600, res = 150)
      print(p)
      dev.off()
      
      plotdata_tmp_df <- plotdata_df %>%
        filter(treatment_length == "1month")
      
      p <- ggscatter(data = plotdata_tmp_df, x = "exp_value", y = "relative_tumor_volume",
                     add = "reg.line",  # Add regressin line
                     label = "model_id",
                     add.params = list(color = "grey", fill = "lightgray", linetype = 2), # Customize reg. line
                     conf.int = T # Add confidence interval
      )
      p <- p + stat_cor(method = "pearson",
                        label.x = min(plotdata_tmp_df$x_plot, na.rm = T),
                        label.y = (max(plotdata_tmp_df$y_plot, na.rm = T) + 0.02))
      p <- p + ggtitle(label = paste0(gene_tmp, " mRNA ~", treatment_tmp, " 1-month response"))
      p <- p + xlab(paste0(gene_tmp, " gene expression at baseline"))
      p <- p + ylab(paste0("Relative tumor volume (", treatment_tmp, "/Control)"))
      p <- p + xlim(c(min(plotdata_tmp_df$x_plot)-0.2, max(plotdata_tmp_df$x_plot)+0.2))
      p <- p + ylim(c(min(plotdata_tmp_df$y_plot)-0.05, max(plotdata_tmp_df$y_plot)+0.05))
      
      file2write <- paste0(dir_tmp_out, gene_tmp, ".", pathwayname_tmp, ".", "1month", ".png")
      png(file2write, width = 600, height = 500, res = 150)
      print(p)
      dev.off()
      
    }
  }
}
