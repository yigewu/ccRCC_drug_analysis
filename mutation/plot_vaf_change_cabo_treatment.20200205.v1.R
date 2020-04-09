# Yige Wu @ WashU 2020 Feb
## plot VAF changes before and after treatment

# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/ccRCC_Drug_shared.R")
library(ggrepel)

# set run id --------------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)

# set output directory ----------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input batch 8 mutation calls --------------------------------------------
mut_batch8_tab <- fread("./Ding_Lab/Projects_Current/PDX-WashU/batch8/somaitcMut/pdx/b8.tumorNormal.sm.non-silent.plus.meta2.filteredPDXmut.tsv", data.table = F)
mut_batch8_tab <- mut_batch8_tab %>%
  filter(grepl(pattern = "RESL", x = Tumor_Sample_Barcode))

# set paired tumor barcodes ----------------------------------------------
treated_tumor_barcode <- "TWDE-WUR-014-RESL5D-4975-Cab-DNA"
control_tumor_barcode <- "TWDE-WUR-014-RESL5D-4971-CT-DNA"

# filter maf for the paired sample --------------------------------------------
maf_tab <- mut_batch8_tab %>%
  filter(Tumor_Sample_Barcode %in% c(treated_tumor_barcode, control_tumor_barcode))

# get VAF for mutations in this pair ----------------------------------------------------------
genes_allmut <- unique(maf_tab$Hugo_Symbol)
vaf_mat <- get_somatic_mutation_vaf_matrix(pair_tab = genes_allmut, maf = maf_tab)
vaf_mat[vaf_mat == ""] <- 0
mut_mat <- get_somatic_mutation_detailed_matrix(pair_tab = genes_allmut, maf = maf_tab)
mut_mat <- mut_mat %>%
  rename(control_mutation = control_tumor_barcode) %>%
  rename(treated_mutation = treated_tumor_barcode) %>%
  mutate(mutation = ifelse(control_mutation == "", treated_mutation, control_mutation))
vaf_mat$mutation <- mut_mat$mutation
file2write <- paste0(dir_out, "Somatic_Mutation_VAF_Changes.Cabo_Treated_vs_Control.", run_id, ".tsv")
write.table(x = vaf_mat, file = file2write, quote = F, sep = "\t", row.names = F)

# get genes overlapping between treated and control -----------------------
genes2plot <- intersect(unique(mut_batch8_tab$Hugo_Symbol[mut_batch8_tab$Tumor_Sample_Barcode == treated_tumor_barcode]),
                        unique(mut_batch8_tab$Hugo_Symbol[mut_batch8_tab$Tumor_Sample_Barcode == control_tumor_barcode]))
genes2plot %>% length()

# plot scatterplot --------------------------------------------------------
plot_data <- vaf_mat %>%
  rename(control_vaf = control_tumor_barcode) %>%
  rename(treated_vaf = treated_tumor_barcode)
## make sure the numbers are in numeric format
plot_data$control_vaf <- as.numeric(as.vector(plot_data$control_vaf))
plot_data$treated_vaf <- as.numeric(as.vector(plot_data$treated_vaf))

plot_data <- plot_data  %>%
  mutate(color = ifelse(treated_vaf > (control_vaf + 0.1), "treated_vaf > control_vaf + 10%", 
                        ifelse(treated_vaf < (control_vaf - 0.1), "treated_vaf < control_vaf - 10%", "|treated_vaf - control_vaf| < 10%")))

## plot mutations in SMGs
### filter the data frame to get the ones to highlight
text_data <- plot_data %>%
  filter((Hugo_Symbol %in% SMGs[["CCRCC"]]))

p <- ggplot()
p <- p + geom_point(data = plot_data, mapping = aes(x = control_vaf, y = treated_vaf))
p <- p + geom_abline(slope = 1, linetype = 2)
p <- p + geom_text_repel(data = text_data, mapping = aes(x = control_vaf, y = treated_vaf, label = Hugo_Symbol, color = color))
p <- p + scale_color_manual(values = c("treated_vaf > control_vaf + 10%" = "red", "treated_vaf < control_vaf - 10%" = "blue", "|treated_vaf - control_vaf| < 10%" = "black"))
p <- p + theme(legend.position = "top")
file2write <- paste0(dir_out, "Somatic_Mutation_VAF_Changes.Cabo_Treated_vs_Control.SMGs.", run_id, ".png")
png(file = file2write, width = 800, height = 800, res = 150)
print(p)
dev.off()

## plot mutations in with VAF Increase
### filter the data frame to get the ones to highlight
text_data <- plot_data %>%
  filter(color == "treated_vaf > control_vaf + 10%")

p <- ggplot()
p <- p + geom_point(data = plot_data, mapping = aes(x = control_vaf, y = treated_vaf))
p <- p + geom_abline(slope = 1, linetype = 2)
p <- p + geom_text_repel(data = text_data, mapping = aes(x = control_vaf, y = treated_vaf, label = Hugo_Symbol, color = color))
p <- p + scale_color_manual(values = c("treated_vaf > control_vaf + 10%" = "red", "treated_vaf < control_vaf - 10%" = "blue", "|treated_vaf - control_vaf| < 10%" = "black"))
p <- p + theme(legend.position = "top")
file2write <- paste0(dir_out, "Somatic_Mutation_VAF_Changes.Cabo_Treated_vs_Control.Treatment_Increase_VAF.", run_id, ".png")
png(file = file2write, width = 800, height = 800, res = 150)
print(p)
dev.off()

## plot mutations in with VAF Decrease
### filter the data frame to get the ones to highlight
text_data <- plot_data %>%
  filter(color == "treated_vaf < control_vaf - 10%")

p <- ggplot()
p <- p + geom_point(data = plot_data, mapping = aes(x = control_vaf, y = treated_vaf))
p <- p + geom_abline(slope = 1, linetype = 2)
p <- p + geom_text_repel(data = text_data, mapping = aes(x = control_vaf, y = treated_vaf, label = Hugo_Symbol, color = color))
p <- p + scale_color_manual(values = c("treated_vaf > control_vaf + 10%" = "red", "treated_vaf < control_vaf - 10%" = "blue", "|treated_vaf - control_vaf| < 10%" = "black"))
p <- p + theme(legend.position = "top")
file2write <- paste0(dir_out, "Somatic_Mutation_VAF_Changes.Cabo_Treated_vs_Control.Treatment_Decrease_VAF.", run_id, ".png")
png(file = file2write, width = 800, height = 800, res = 150)
print(p)
dev.off()


