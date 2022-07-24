# Yige Wu @WashU Mar 2022

#  set up libraries and output directory -----------------------------------
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

# input dependencies ------------------------------------------------------
## input protein data
exp_df <- fread(data.table = F, input = "../ccRCC_snRNA/Resources/Bulk_Processed_Data/Protein/6_CPTAC3_CCRCC_Whole_abundance_gene_protNorm=2_CB.tsv")
## input bulk meta data
metadata_bulk_df <- fread("../ccRCC_snRNA/Resources/Bulk_Processed_Data/Case_ID/CPTAC_ccRCC_discovery_caseID_v1.0.tsv")

# preprocess ------------------------------------------------------
metadata_filtered_df <- metadata_bulk_df %>%
  filter(Histologic_Type == "Clear cell renal cell carcinoma") %>%
  filter(!is.na(Specimen.Label.normal))
# > nrow(metadata_filtered_df)
# [1] 81
## center data
sampleids_tumor <- metadata_filtered_df$Specimen.Label.tumor
sampleids_nat <- metadata_filtered_df$Specimen.Label.normal

pos <- position_jitter(width = 0.2, seed = 3)

# for (gene_test in c("TGFBI", "LRP1", "COL7A1", "TAGLN", "SCG2", "THBS1", "TFPI2")) {
# for (gene_test in c("CDK4", "TNPO2", "RPS6KA5", "PRMT5", "NUP98")) {
# for (gene_test in c("TK1", "CDK1", "KPNA2", "MCM2", "PRDX4", "PCNA", "TBRG4", "MCM3", "RFC3", "MCM7")) {
# for (gene_test in c("TIMM8B", "GOT2", "HCCS", "UQCRQ", "SLC25A20", "NDUFB1", "COX7C", "TOMM22", "SLC25A5", "UQCRB", "UQCR10", "COX4I1", "COX17", "MRPL11", "DLD")) {
# for (gene_test in c("EIF3A", "MRPS9", "PTCD3", "RPL13A", "EIF4G1", "TSFM", "MRPL41", "RPS5", "MRPL14", "MRPL16", "MRPL9", "MRPL3", "MRPL19", "MRPS31", "MRPL11", "MRPL43", "MRPL1", "MRPS10")) {
for (gene_test in c("CHMP6", "IGF2BP3", "DNAJC7", "PYCR1", "MRM3", "ERO1B", "CLIC6", "COBLL1")) {
  # for (gene_test in c("ICAM1", "LRP1", "FTH1", "ICAM1", "ALDH1L1", "CDK4", "BLVRB", "ADK", "MAT1A")) {
  ## filter specific protein data
  exp_tumor_df <- exp_df[exp_df$Index == gene_test, sampleids_tumor]
  testdata_tumor_df <- data.frame(CASE_ID = metadata_filtered_df$CASE_ID, Expression = unlist(exp_tumor_df), Group = "Tumor")
  
  exp_nat_df <- exp_df[exp_df$Index == gene_test, sampleids_nat]
  testdata_nat_df <- data.frame(CASE_ID = metadata_filtered_df$CASE_ID, Expression = unlist(exp_nat_df), Group = "NAT")
  
  plotdata_df <- rbind(testdata_tumor_df, testdata_nat_df)
  plotdata_df <- plotdata_df %>%
    filter(!is.na(Expression))
  
  p <- ggplot(data = plotdata_df, mapping = aes(x = Group, y = Expression))
  p <- p + geom_violin(mapping = aes(fill = Group))
  p = p + geom_boxplot(width=.1, outlier.shape = 23, outlier.fill = "black")
  # p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.8, size = 1.5)
  p = p + stat_compare_means(data = plotdata_df, 
                             mapping = aes(x = Group, y = Expression),
                             # symnum.args = symnum.args,
                             hide.ns = F, method = "wilcox.test", size = 7, 
                             label.y = max(plotdata_df$Expression) + 0.1*(max(plotdata_df$Expression)-min(plotdata_df$Expression)))
  p <- p + ylim(c(min(plotdata_df$Expression), max(plotdata_df$Expression) + 0.2*(max(plotdata_df$Expression)-min(plotdata_df$Expression))))
  p <- p + scale_fill_manual(values = c("NAT" = "#4DAF4A", "Tumor" = "#E41A1C"))
  p <- p + theme_classic(base_size = 18)
  p <- p + theme(legend.position = "none")
  p <- p + ylab(label = paste0("Bulk protein level\n(CPTAC ccRCC cohort)")) + ggtitle(label = paste0(gene_test, " protein level"))
  p <- p + theme(axis.title.x = element_blank(), 
                 axis.title.y = element_text(size = 18), 
                 axis.text = element_text(color = "black", size = 18))
  
  file2write <- paste0(dir_out, gene_test, ".pdf")
  pdf(file2write, width = 4.5, height = 4, useDingbats = F)
  print(p)
  dev.off()
}
