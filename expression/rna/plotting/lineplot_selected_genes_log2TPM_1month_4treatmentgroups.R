# Yige Wu @ WashU 2023 Jun

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("OmnipathR")
packages = c(
  "plyr",
  "stringr",
  "reshape2",
  "data.table",
  "dplyr"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

# input dependencies ------------------------------------------------------
## input the protein data
exp_df <- fread("./Data_Freeze/v1.dataFreeze.washU_rcc/3.geneExp/v3.20210116/datafreeze.v3.kallisto.geneExp.protein_coding.tsv", data.table = F)
## input the sample info
sampleinfo_df <- readxl::read_excel(path = "./Data_Freeze/v1.dataFreeze.washU_rcc/0.sample_info/v3.20210116/RCC_PDX_Samples.20210115.v2.xlsx")

# set parameters ----------------------------------------------------------
genes_filter <- c("CCND3", "MYBL2", "FOXM1", "CDK2", "BCL2L12","BCL2L13")
genes_filter <- c("CCND3", "MYBL2", "FOXM1", "BCL2L12")
colnames_id <- colnames(exp_df)[!(colnames(exp_df) %in% sampleinfo_df$Analysis_ID)]

# plot --------------------------------------------------------------------
plot_data_wide_df <- exp_df %>%
  filter(Name %in% genes_filter)
plot_data_long_df <- melt(data = plot_data_wide_df, id.vars = colnames_id)
plot_data_long_df <- as.data.frame(plot_data_long_df)
plot_data_long_df <- merge(x = plot_data_long_df, y = sampleinfo_df, by.x = "variable", by.y = "Analysis_ID", all.x = T)

plot_data_long_df <- plot_data_long_df %>%
  filter(ShortTag %in% c("Control", "Treated.Cabo+Sap")) %>%
  filter(Treatment.Month == 1) %>%
  filter(!(PairTag %in% c("RESL5D.27", "RESL5E.28", "RESL5E.30", "RESL10F.38"))) %>%
  mutate(id_model_length = paste0(ifelse(is.na(PairTag), ModelID, PairTag), "_", Treatment.Month))
plot_data_long_df$Treatment = mapvalues(plot_data_long_df$ShortTag, 
                                        from = c("Control", "Treated.Cabo", "Treated.Sap", "Treated.Cabo+Sap"),
                                        to = c("Control", "Cabozantinib", "Sapanisertib", "Cabozantinib+\nSapanisertib"))
ids_process = unique(plot_data_long_df$id_model_length[plot_data_long_df$Treatment == "Cabozantinib+\nSapanisertib"])
plot_data_long_df = plot_data_long_df %>%
  filter(id_model_length %in% ids_process) %>%
  arrange(Name, id_model_length, desc(Treatment))

plot_data_long_df$Treatment = factor(x = plot_data_long_df$Treatment, levels = c("Control", "Cabozantinib", "Sapanisertib", "Cabozantinib+\nSapanisertib"))

plot_data_long_df$Model_id <- factor(x =   plot_data_long_df$ModelID, levels = c("RESL5", "RESL10", "RESL11", "RESL3", "RESL12", "RESL4"))
plot_data_long_df$y_plot = log2(plot_data_long_df$value + 1)
# plot_data_long_df$y_plot <- (plot_data_long_df$value - min(plot_data_long_df$value))/(max(plot_data_long_df$value) - min(plot_data_long_df$value))

stat.test <- plot_data_long_df %>%
  group_by(Name) %>%
  t_test(y_plot ~ Treatment, paired = T, 
         comparisons = list(c("Cabozantinib+\nSapanisertib", "Control"))) %>%
  adjust_pvalue() %>%
  add_significance("p")
# stat.test <- stat.test %>% add_xy_position(fun = "mean_se", x = "Treatment")
stat.test <- stat.test %>% add_xy_position(fun = "mean_se", x = "Treatment")
genes_sorted = stat.test %>%
  arrange(p)
plot_data_long_df$Name = factor(x = plot_data_long_df$Name, levels = genes_sorted$Name)
plot_data_long_df$Name = factor(x = plot_data_long_df$Name, levels = genes_filter)

## plot
p <- ggline(data = plot_data_long_df, 
            x = "Treatment", y = "y_plot", 
            color = "Model_id", 
            shape = "Treatment", size = 0.6, add.params = list(alpha = 0.6),
            position = position_dodge(),
            facet.by = "Name", nrow = 1)
# p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01)
p <- p + labs(y = paste0(" gene expression\n(log2(TPM+1))"))
p <- p + scale_color_discrete(name = "") + 
  guides(color=guide_legend(ncol=2, byrow=FALSE)) 
  # guides(color=guide_legend(nrow=2, byrow=TRUE)) 
p <- p + scale_shape_discrete(name = "")
p <- p + theme_classic(base_size = 18)
p <- p + theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.x = element_blank(),
               # legend.box="vertical", 
               legend.margin = margin(),
               legend.position = "right")

# write output ------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

## write output
file2write <- paste0(dir_out, "lineplot", ".pdf")
pdf(file2write, width = 8, height = 3, useDingbats = F)
print(p)
dev.off()
# file2write <- paste0(dir_out, gene_plot, ".png")
# png(file2write, width = 1500, height = 800, res = 150)
# print(p)
# dev.off()
