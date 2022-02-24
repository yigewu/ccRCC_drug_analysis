# Yige Wu @ WashU 2022 Jan

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
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
pathwayscores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/calculate_pathway_score_treatedvscontrol_bydiff_log2intensity/20220112.v3/Pathway_scores_bysample.20220112.v3.tsv")
## input the relative tumor volume
tumorvolume_df <- readxl::read_excel(path = "./Resources/Treatment_Lists/Treatment_Summary/Cabo_Sapa_TumorReduction.xlsx", sheet = "1month_2month_divide")

# plot by pathway ---------------------------------------------------------
pathwayname_tmp <- "Cell_cycle"
treatment_tmp <- "Sap"
for (treatment_tmp in c("Cabo+ Sap")) {
# for (treatment_tmp in c("Cabo", "Sap", "Cabo+ Sap")) {
  for (pathwayname_tmp in c("Cell_cycle", "Apoptosis", "RTK", "PI3K_Akt", "Ras_MAPK", "EMT", "DDR")) {
    sampleids_tmp <- colnames(pathwayscores_df)
    sampleids_tmp <- sampleids_tmp[grepl(x = sampleids_tmp, pattern = paste0("_", ifelse(treatment_tmp == "Cabo+ Sap", "Cabo\\+ Sap", treatment_tmp), "_"))]
    pathwayscores_tmp <- unlist(pathwayscores_df[pathwayscores_df$Pathway_name == pathwayname_tmp, sampleids_tmp])
    tumorvolumes_tmp <- mapvalues(x = sampleids_tmp, from = paste0(tumorvolume_df$Model, "_", treatment_tmp, "_", tumorvolume_df$Treatment_length), to = unlist(tumorvolume_df[,paste0(treatment_tmp, " vs. Control")]))
    tumorvolumes_tmp <- as.numeric(tumorvolumes_tmp)
    
    ## make plot data
    plotdata_df <- data.frame(sample_id = sampleids_tmp,
                              model_id = str_split_fixed(string = sampleids_tmp, pattern = "_", n = 3)[,1],
                              treatment_length = str_split_fixed(string = sampleids_tmp, pattern = "_", n = 3)[,3],
                              pathway_score = pathwayscores_tmp, 
                              relative_tumor_volume = tumorvolumes_tmp)
    
    ## plot
    p <- ggscatter(data = plotdata_df, x = "pathway_score", y = "relative_tumor_volume", shape = "treatment_length",
                   add = "reg.line",  # Add regressin line
                   # add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                   conf.int = TRUE # Add confidence interval
    )
    p <- p + stat_cor(method = "pearson",
                      label.x = quantile(x = plotdata_df$pathway_score, probs = 0.001, na.rm = T))
    p <- p + ggtitle(label = paste0(pathwayname_tmp, "~", treatment_tmp, " treatment"))
    p <- p + xlab(paste0(pathwayname_tmp, " score (", treatment_tmp, " vs. Control)"))
    p <- p + ylab(paste0("Relative tumor volume (", treatment_tmp, "/Control)"))
    file2write <- paste0(dir_out, pathwayname_tmp, "_", treatment_tmp, ".png")
    png(file2write, width = 600, height = 600, res = 150)
    print(p)
    dev.off()
    
    p <- ggscatter(data = plotdata_df, x = "pathway_score", y = "relative_tumor_volume", shape = "treatment_length", color = "treatment_length",
                   add = "reg.line",  # Add regressin line
                   # add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                   conf.int = TRUE # Add confidence interval
    )
    p <- p + stat_cor(aes(color = treatment_length),
                      method = "pearson",
                      label.x = quantile(x = plotdata_df$pathway_score, probs = 0.001, na.rm = T))
    p <- p + ggtitle(label = paste0(pathwayname_tmp, "~", treatment_tmp, " treatment"))
    p <- p + xlab(paste0(pathwayname_tmp, " score (", treatment_tmp, " vs. Control)"))
    p <- p + ylab(paste0("Relative tumor volume (", treatment_tmp, "/Control)"))
    file2write <- paste0(dir_out, pathwayname_tmp, "_", treatment_tmp, ".bytreatmentlength.png")
    png(file2write, width = 600, height = 600, res = 150)
    print(p)
    dev.off()
    
  }
}
