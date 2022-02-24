# Yige Wu @ WashU 2022 Jan

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
source("./ccRCC_drug_analysis/plotting.R")
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
tumorvolume_df <- readxl::read_excel(path = "./Resources/Treatment_Lists/Treatment_Summary/Cabo_Sapa_TumorReduction.xlsx", sheet = "1-month_divide")

# plot by pathway ---------------------------------------------------------
pathwayname_tmp <- "Cell_cycle"
treatment_tmp <- "Sap"
sampleids_tmp <- colnames(pathwayscores_df)
sampleids_tmp <- sampleids_tmp[grepl(x = sampleids_tmp, pattern = paste0("_", treatment_tmp, "_")) & grepl(x = sampleids_tmp, pattern = "1month")]
modelids_tmp <- str_split_fixed(string = sampleids_tmp, pattern = "_", n = 3)[,1]
  
pathwayscores_tmp <- unlist(pathwayscores_df[pathwayscores_df$Pathway_name == pathwayname_tmp, sampleids_tmp])
tumorvolumes_tmp <- mapvalues(x = modelids_tmp, from = tumorvolume_df$Model, to = unlist(tumorvolume_df[,paste0(treatment_tmp, " vs. Control")]))
tumorvolumes_tmp <- as.numeric(tumorvolumes_tmp)

## make plot data
plotdata_df <- data.frame(model_id = modelids_tmp, pathway_score = pathwayscores_tmp, relative_tumor_volume = tumorvolumes_tmp)

## plot
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = pathway_score, y = relative_tumor_volume))
p
