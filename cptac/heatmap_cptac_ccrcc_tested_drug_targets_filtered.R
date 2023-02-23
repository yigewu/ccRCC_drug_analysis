# Yige Wu @ WashU 2019 Nov
## plot a heatmap with proteomics data from the discovery set data freeze


# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2",
  "ComplexHeatmap",
  "circlize"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
source("./ccRCC_drug_analysis//functions.R")
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input data --------------------------------------------------------------
## input cobined protein and phosphoprotein data
exp_df <-fread(data.table = F, input = "./Resources/Analysis_Results/cptac/unite_cptac_ccRCC_protein_w_phospho_data/20221216.v1/Protein_Phosphorylation_combined.refnormalized.20221216.v1.tsv")
## input bulk meta data
bulk_meta_tab <- fread("../../CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/cptac-metadata.csv")
## input HIF1A targets
hif_targets_df <- fread(data.table = F, input = "../ccRCC_snRNA/Resources/Analysis_Results/dependencies/write_hif_targets/20201027.v1/HIF_Target_Genes.20201027.v1.tsv")
# ## input tumor vs normal data
protein_test_df <- fread(data.table = F, input = "../ccRCC_snRNA/Resources/Analysis_Results/bulk/expression/protein/compare_bulk_protein_tumor_vs_normal/20221219.v1/Bulk_Protein_Tumor_vs_Normal.Wilcox.20221219.v1.tsv")
phospho_test_df <- fread(data.table = F, input = "../ccRCC_snRNA/Resources/Analysis_Results/bulk/expression/phosphorylation/compare_bulk_phosphorylation_tumor_vs_matched_NATs/20221219.v1/Bulk_Phosphosite_Tumor_vs_NAT.Matched.Wilcox.20221219.v1.tsv")

# input parameters --------------------------------------------------------
## set plotting parameters
num_nonna <- 0
row_fontsize <- 10
## set header names
header_colnames <- c("Protein_id", "Gene", "Phosphosite", "Index", "ReferenceIntensity")
aliquotids_colnames <- colnames(exp_df)[!(colnames(exp_df) %in% header_colnames)]


# prepare phosphosites and proteins to plot -------------------------------
## mTOR inhibitor: https://www.researchgate.net/figure/The-mammalian-target-of-rapamycin-complex-1-mTORC1-pathway-and-translation-initiation_fig1_309236845
## reference: https://www.sciencedirect.com/science/article/pii/S0092867417301824?via%3Dihub#bib19
## CDK inhibitor: https://pubchem.ncbi.nlm.nih.gov/compound/Abemaciclib#:~:text=Abemaciclib%20is%20an%20orally%20available,protein%20phosphorylation%20in%20early%20G1.
## IAP inhibitor: https://www.frontiersin.org/articles/10.3389/fonc.2020.532292/full#:~:text=Birinapant%20is%20a%20SMAC%20mimetic,)%20and%20XIAP%20(10).
hif1_targets_vec <- hif_targets_df$target_genesymbol[hif_targets_df$source_genesymbol == "HIF1A" & !is.na(hif_targets_df$target_genefunction)]
proteinids_plot_df <- data.frame(Gene = c("MET", "MET", "AXL", "FLT1", "KDR", "RET", "KIT", "TIE2",
                                          "AKT2", "AKT2", "RPS6", rep("EIF4EBP1", 5),
                                          "BIRC3", "BIRC2", "XIAP",
                                          "HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC5", "HDAC6", "HDAC7", "HDAC8", "HDAC9", "HDAC10", "HDAC11", 
                                          "CDK4", "CDK6",
                                          "MAP2K1", "MAP2K1", "MAP2K1", "MAPK1", "MAPK3",
                                          "HIF1A", hif1_targets_vec,
                                          "AGTR1", "AGT"),
                                 Phosphosite = c("Y1234", rep("Protein", 7),
                                                 "T309", "S474", "S235S236", "T37", "T46", "T70", "T37T46", "S65",
                                                 rep("Protein", 3),
                                                 rep("Protein", 11),
                                                 rep("Protein", 2),
                                                 "S218", "S222", "S218S222", "T185Y187", "T202Y204",
                                                 "Protein", rep(x = "Protein", length(hif1_targets_vec)),
                                                 rep("Protein", 2)),
                                 Pathway = c(rep("RTKs", 8),
                                             rep("mTOR\nsignaling", 8),
                                             rep("IAPs", 3),
                                             rep("HDACs", 11),
                                             rep("CDKs", 2),
                                             rep("MEK/ERK", 5),
                                             rep("HIF\nsignaling", length(hif1_targets_vec)+1),
                                             rep("Angiotensin\nreceptor", 2)))
proteinids_plot_df <- proteinids_plot_df %>%
  mutate(Protein_id = paste0(Gene, "_", Phosphosite))

# process aliquot ids --------------------------------------------
### get the ids for the normal samples
normal_bulk_aliquot_ids2plot <- bulk_meta_tab$Specimen.Label[bulk_meta_tab$Set.A == "yes" & bulk_meta_tab$Type == "Normal"]
normal_bulk_aliquot_ids2plot
case_ids2plot <- mapvalues(x = normal_bulk_aliquot_ids2plot, from = bulk_meta_tab$Specimen.Label, to = bulk_meta_tab$Case.ID)
case_ids2plot
tumor_bulk_aliquot_ids2plot <- mapvalues(x = case_ids2plot, from = bulk_meta_tab$Case.ID[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"], to = as.vector(bulk_meta_tab$Specimen.Label[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"]))
tumor_bulk_aliquot_ids2plot

# filter and format  expression data -------------------------------------
## filter the phosphosite
plot_data_df <- exp_df %>%
  filter(Protein_id %in% proteinids_plot_df$Protein_id)
## remove duplicates
plot_data_df <- plot_data_df[!duplicated(plot_data_df$Protein_id),]
## take out the reference intensity
plot_data_mat <- as.matrix(plot_data_df[,aliquotids_colnames])
rownames(plot_data_mat) <- plot_data_df$Protein_id
## add information for each protein
number_values_vec <- rowSums(!is.na(plot_data_mat))
proteinids_plot_df$number_value <- number_values_vec[proteinids_plot_df$Protein_id]
# proteins_sigup <- protein_test_df$symbol[!is.na(protein_test_df$pvaladj) & protein_test_df$pvaladj < 0.05 & protein_test_df$logFC.ToverN > 0]; proteins_sigup <- paste0(proteins_sigup, "_Protein")
proteins_sigup <- protein_test_df$gene_symbol[!is.na(protein_test_df$fdr) & protein_test_df$fdr < 0.05 & protein_test_df$meddiff_exp > 0]; proteins_sigup <- paste0(proteins_sigup, "_Protein")
phospho_test_df <- phospho_test_df %>% mutate(Protein_id = paste0(SUBSTRATE, "_", SUB_MOD_RSD))
phosphosites_sigup <- phospho_test_df$Protein_id[!is.na(phospho_test_df$fdr) & phospho_test_df$fdr < 0.05 & phospho_test_df$meddiff_exp > 0]
proteinids_plot_df <- proteinids_plot_df %>%
  mutate(Upregulated_in_tumor = ifelse(Protein_id %in% c(proteins_sigup, phosphosites_sigup), "Yes", "No"))
## filter proteins
proteinids_plot_filtered_df <- proteinids_plot_df %>%
  filter(!is.na(number_value)) %>%
  filter((number_value >= num_nonna & Upregulated_in_tumor == "Yes") | Protein_id %in% c("HIF1A_Protein", "AGTR1_Protein", "MET_Y1234"))
plot_data_mat2 <- plot_data_mat[proteinids_plot_filtered_df$Protein_id,]
## make row ids
rowids_plot <- rownames(plot_data_mat2)
## clean up
rm(exp_df)

# make colors -------------------------------------------------------------
col_protein = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# make column annotation --------------------------------------------------
# ca = HeatmapAnnotation(Sample_Type = ifelse(colnames(plot_data_mat2) %in% tumor_bulk_aliquot_ids2plot, "Tumor", "Normal"),
#                        col = list(Sample_Type = c("Tumor" = "#E41A1C", "Normal" = "#4DAF4A")))

# make column annotation --------------------------------------------------
ra = rowAnnotation(Upregulated_in_Tumor = mapvalues(x = rowids_plot, from = proteinids_plot_df$Protein_id, to = as.vector(proteinids_plot_df$Upregulated_in_tumor)),
                       col = list(Upregulated_in_Tumor = c("Yes" = "red", "No" = "black")), annotation_name_side = "top", annotation_name_rot = 0)

# make column split -------------------------------------------------------
col_split_vec <- ifelse(colnames(plot_data_mat2) %in% tumor_bulk_aliquot_ids2plot, "Tumor", "Normal")

# make row split ----------------------------------------------------------
row_split_vec <- mapvalues(x = rowids_plot, from = proteinids_plot_df$Protein_id, to = as.vector(proteinids_plot_df$Pathway))
row_split_factor <- factor(x = row_split_vec, levels = unique(proteinids_plot_df$Pathway))

# plot --------------------------------------------------------------------
p <- Heatmap(matrix = plot_data_mat2,
             name = "log2Intensity\n(Sample/Reference)", 
             col = col_protein,
             ## column
             # top_annotation = ca, 
             column_split = col_split_vec,
             cluster_columns = T, show_column_dend = F, cluster_column_slices = T,
             show_column_names = F,
             ## row
             right_annotation = ra,
             row_title_gp = gpar(fontsize = 12), row_title_rot = 0,
             row_names_gp = gpar(fontsize = row_fontsize), row_split = row_split_factor,
             cluster_rows = F)
p
## save heatmap to PNG
file2write <- paste0(dir_out, "CPTAC.", run_id, ".pdf")
pdf(file2write, width = 10, height = 8, useDingbats = F)
print(p)
dev.off()


