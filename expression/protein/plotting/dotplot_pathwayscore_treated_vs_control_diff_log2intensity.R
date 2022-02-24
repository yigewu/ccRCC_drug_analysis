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
version_tmp <- 3
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input pathway score changes for treated vs. control
# pathwayscores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/calculate_pathway_score_treatedvscontrol_bydiff_log2intensity/20220111.v1/Pathway_scores_bysample.20220111.v1.tsv")
# pathwayscores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/calculate_pathway_score_treatedvscontrol_bydiff_log2intensity/20220112.v1/Pathway_scores_bysample.20220112.v1.tsv")
pathwayscores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/other/calculate_pathway_score_treatedvscontrol_bydiff_log2intensity/20220112.v3/Pathway_scores_bysample.20220112.v3.tsv")

# set plotting parameters -------------------------------------------------
treatment_tmp <- "Cabo"
sample_ids <- colnames(pathwayscores_df)
sample_ids <- sample_ids[grepl(x = sample_ids, pattern = "RESL")]
# sample_ids <- sample_ids[grepl(x = sample_ids, pattern = "RESL") & grepl(x = sample_ids, pattern = paste0("_", treatment_tmp, "_"))]
sample_ids <- sample_ids[grepl(x = sample_ids, pattern = "RESL5|RESL4|RESL10")]
sample_ids
colors_bymodel <- RColorBrewer::brewer.pal(n = 7, name = "Set1")[-6]
names(colors_bymodel) <- c("RESL5", "RESL10", "RESL12", "RESL4", "RESL11", "RESL3")

# make plot data ----------------------------------------------------------
plotdata_df <- melt(pathwayscores_df[, c("Pathway_name", sample_ids)])
plotdata_df <- plotdata_df %>%
  mutate(model_id = str_split_fixed(string = variable, pattern = "_", n = 3)[,1]) %>%
  mutate(treatment = str_split_fixed(string = variable, pattern = "_", n = 3)[,2]) %>%
  mutate(treatment_length = str_split_fixed(string = variable, pattern = "_", n = 3)[,3]) %>%
  filter(!is.na(value)) %>%
  mutate(y_plot = ifelse(value > 1, 1,
                         ifelse(value < -1, -1, value)))
plotdata_df$treatment <- factor(x = plotdata_df$treatment, levels = c("Cabo", "Sap", "Cabo+ Sap"))
plotdata_df$Pathway_name <- factor(x = plotdata_df$Pathway_name, levels = c("Cell_cycle", "RTK", "PI3K_Akt", "Ras_MAPK",  "Apoptosis", "EMT", "DDR"))

# plot --------------------------------------------------------------------
p <- ggplot()
# p <- p + geom_violin(data = plotdata_df, mapping = aes(x = Pathway_name, y = y_plot), color = NA, fill = "grey50", alpha = 0.5)
p <- p + geom_dotplot(data = plotdata_df, mapping = aes(x = Pathway_name, y = y_plot, fill = model_id), 
                      binaxis = "y", stackdir = "center", alpha = 0.7, position=position_dodge(0.3))
# p <- p + scale_shape_manual(values = c("1month" = 20, "2month" = 19))
p <- p + scale_fill_manual(values = colors_bymodel)
p <- p + facet_grid(treatment_length~treatment)
p <- p + theme_bw()
p <- p + ylab("Pathway activity change by Cabo/Sap inhibition")
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p <- p + theme(axis.title.x = element_blank())
# p
file2write <- paste0(dir_out, "Pathway_score_treated_vs_control.", run_id, ".pdf")
pdf(file2write, width = 7, height = 4, useDingbats = F)
print(p)
dev.off()
file2write <- paste0(dir_out, "Pathway_score_treated_vs_control.", run_id, ".png")
png(file2write, width = 2000, height = 800, res = 150)
print(p)
dev.off()
