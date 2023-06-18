# Yige Wu @WashU Jun 2020

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
  "readxl",
  "ggplot2",
  "ggpubr"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg_name_tmp)
  }
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_drug_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input formatted measurement data --------------------------------------------------
measurement_df <- fread(data.table = F, input = "./Resources/Analysis_Results/treatment_data/extract_treated_pdx_measurement/20200616.v1/RCC_PDX.Measurement.Long.20200616.v1.tsv")

# extract body weight and get starting and 1-month body weight ------------
measurement_df <- measurement_df %>%
  mutate(ModelID = str_split_fixed(string = Tab, pattern = "_", n = 3)[,1]) %>%
  mutate(Batch = str_split_fixed(string = Tab, pattern = "_", n = 3)[,2]) %>%
  mutate(PDXLineID = paste0(ModelID, "_", Batch))
unique(measurement_df$PDXLineID)
## filter PDX line by those with complete cohort of 4 Cabo+Sap single+combo groups
measurement_filtered <- measurement_df %>%
  filter(PDXLineID %in% c("RESL5_B2", "RESL10_B1", "RESL12B_B1", "RESL4_B1", "RESL11D_B1", "RESL3_B1"))
unique(measurement_filtered$Group)
## filter by tratment group
measurement_filtered <- measurement_filtered %>%
  filter(Group %in% c("CT", "Cab", "Sap", "Cab+Sap"))
View(measurement_filtered %>%
       select(PDXLineID, Group) %>%
       unique() %>%
       table())
## get body weight only
bw_df <- measurement_filtered %>%
  filter(Measurement_Type == "B.W./g") %>%
  mutate(Bodyweight = as.numeric(value))
## get body weight at treatment start
bw_treatmentstart_df <- bw_df %>%
  filter(Date == Date.TreatmentStart) %>%
  mutate(Bodyweight.TreatmentStart = Bodyweight) %>%
  select(PDXLineID, ID.Xiaolu, Group, Bodyweight.TreatmentStart)
bw_treatmentstart_df$Bodyweight.TreatmentStart
## get body weight at closest to 1-month
bw_treatment1month_df <- bw_df %>%
  group_by(PDXLineID, ID.Xiaolu, Group) %>%
  summarize(Days.ClosestTo1Month = Treatment.Days[which.min(abs(Treatment.Days - 30))], Bodyweight.ClosestTo1Month = Bodyweight[which.min(abs(Treatment.Days - 30))])
bw_treatment1month_df <- bw_treatment1month_df %>%
  filter(Days.ClosestTo1Month >= 23)
bw_treatment1month_df$Bodyweight.ClosestTo1Month

# make plot data -----------------------------------------------------------
plot_data_df <- merge(bw_treatmentstart_df,
                      bw_treatment1month_df,
                      by = c("PDXLineID", "ID.Xiaolu", "Group"))
plot_data_df <- plot_data_df %>%
  mutate(BodyweightChange = (Bodyweight.ClosestTo1Month - Bodyweight.TreatmentStart)/Bodyweight.TreatmentStart) %>%
  mutate(plot_x = Group) %>%
  mutate(plot_y = BodyweightChange*100) %>%
  mutate(model = str_split_fixed(PDXLineID, "_", 2)[,1])
plot_data_df$plot_x <- factor(x = plot_data_df$plot_x, levels = c("CT", "Cab", "Sap", "Cab+Sap"))

# make colors -------------------------------------------------------------
RColorBrewer::display.brewer.all()
RColorBrewer::display.brewer.pal(name = "Set1", n = 7)
colors_group <- c("grey50", RColorBrewer::brewer.pal(name = "Set1", n = 6)[c(1,3,6)])
names(colors_group) <- c("CT", "Cab", "Sap", "Cab+Sap")

# make box plot without annotating model -----------------------------------------------------------
p <- ggplot(data = plot_data_df, mapping = aes(x = plot_x, y = plot_y))
p <- p + geom_violin(mapping = aes(fill = plot_x))
p = p + geom_boxplot(width=.1, outlier.shape = 23, outlier.fill = "black")

p = p + stat_compare_means(comparisons = list(c("Cab", "CT"), c("Sap", "CT"), c("Cab+Sap", "CT")),
                          # symnum.args = symnum.args,
                           hide.ns = F, method = "wilcox.test")
p <- p + geom_hline(yintercept = 0, linetype = 2, alpha = 0.8)
p <- p + scale_fill_manual(values = colors_group)
p <- p + theme_bw()
p <- p + ylab(label = "% Body Weight Change (1-Month)")
p <- p + theme(legend.position = "none",
               axis.title.x = element_blank(),
               axis.text = element_text(size = 15, color = "black"),
               axis.title.y = element_text(size = 15))

file2write <- paste0(dir_out, "BodyWeightChange.", "1MonthTreatment.", "pdf")
pdf(file2write, width = 4, height = 3.5, useDingbats = F)
print(p)
dev.off()
file2write <- paste0(dir_out, "BodyWeightChange.", "1MonthTreatment.", "png")
png(filename = file2write, width = 1000, height = 800, res = 150)
print(p)
dev.off()

# make box plot with annotating model -----------------------------------------------------------
library(ggrepel)
library(RColorBrewer)
p <- ggplot(data = plot_data_df, mapping = aes(x = plot_x, y = plot_y))
# p <- p + geom_violin(trim = F)
p = p + geom_boxplot(width = 0.5)
p = p + geom_dotplot(binaxis = 'y', stackdir = "center", stackgroups = T,
                     binpositions="all",
                     mapping = aes(fill = model))
# p = p + stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1))
p = p + stat_compare_means(comparisons = list(c("Cab", "CT"), c("Sap", "CT"), c("Cab+Sap", "CT")),
                           # symnum.args = symnum.args,
                           hide.ns = F, method = "wilcox.test")
p <- p + geom_hline(yintercept = 0, linetype = 2, alpha = 0.8)
p <- p + scale_fill_manual(values = RColorBrewer::brewer.pal(n = 6, name = "Set1"))
p <- p + theme_bw()
p <- p + ylab(label = "% Body Weight Change (1-Month)")
p <- p + theme(axis.title.x = element_blank(),
               axis.text = element_text(size = 15, color = "black"),
               axis.title.y = element_text(size = 15))
p
file2write <- paste0(dir_out, "BodyWeightChange.", "1MonthTreatment.markmodel", ".pdf")
pdf(file2write, width = 4, height = 3.5, useDingbats = F)
print(p)
dev.off()
file2write <- paste0(dir_out, "BodyWeightChange.", "1MonthTreatment.markmodel", ".png")
png(filename = file2write, width = 700, height = 500, res = 150)
print(p)
dev.off()



