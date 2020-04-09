# Yige Wu @ WashU 2019 Nov
## plot a heatmap with proteomics data from the discovery set data freeze


# set up libraries and output directory -----------------------------------
## set up working directory and source functions and load libraries
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/ccRCC_Drug_shared.R")
## set run id 
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input dependencies --------------------------------------------
## set plotting parameters
num_nonna <- 40
row_fontsize <- 9
## set genes to plot
## reference: https://www.researchgate.net/figure/The-mammalian-target-of-rapamycin-complex-1-mTORC1-pathway-and-translation-initiation_fig1_309236845
## reference: https://www.sciencedirect.com/science/article/pii/S0092867417301824?via%3Dihub#bib19
phosphosites <- data.frame(gene = c("RPS6KB1", rep("RPS6", 2), "EIF4B",
                                    rep("EIF4EBP1", 4), "EIF4E"),
                           phosphosite = c("T389", "S236", "S235", "S422",
                                           "T37", "T46","T70", "S65", "S209"))

## input phospho data
pho_df <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/phosphoproteome/6_CPTAC3_CCRCC_Phospho_abundance_phosphosite_protNorm=2_CB.tsv", data.table = F)
## input bulk meta data
bulk_meta_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/cptac-metadata.csv")
## get the bulk aliquot IDs
### get the ids for the normal samples
normal_bulk_aliquot_ids2plot <- bulk_meta_tab$Specimen.Label[bulk_meta_tab$Set.A == "yes" & bulk_meta_tab$Type == "Normal"]
normal_bulk_aliquot_ids2plot
case_ids2plot <- mapvalues(x = normal_bulk_aliquot_ids2plot, from = bulk_meta_tab$Specimen.Label, to = bulk_meta_tab$Case.ID)
case_ids2plot
tumor_bulk_aliquot_ids2plot <- mapvalues(x = case_ids2plot, from = bulk_meta_tab$Case.ID[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"], to = as.vector(bulk_meta_tab$Specimen.Label[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"]))
tumor_bulk_aliquot_ids2plot

# filter and format  expression data -------------------------------------
## make the matrix to plot the heatmap
pho_df %>% head()
header_col_names_keep <- c("Gene", "phosphosite", "Index", "ReferenceIntensity")
plot_data_df <- pho_df %>%
  mutate(phosphosite = str_split_fixed(string = Index, pattern = "_", n = 7)[,7]) %>%
  filter(phosphosite != "")  %>%
  select(header_col_names_keep, tumor_bulk_aliquot_ids2plot, normal_bulk_aliquot_ids2plot)
## filter the phosphosite
row_numbers_filter <- sapply(1:nrow(phosphosites), function(i, phosphosites_df, phosphosite_data_df) {
  gene_tmp <- phosphosites_df$gene[i]
  phosphosite_tmp <- phosphosites_df$phosphosite[i]
  row_numbers <- which(phosphosite_data_df$Gene == gene_tmp & grepl(x = phosphosite_data_df$phosphosite, pattern = phosphosite_tmp))
  return(row_numbers)
}, phosphosites_df = phosphosites, phosphosite_data_df = plot_data_df)
row_numbers_filter <- unlist(row_numbers_filter)
row_numbers_filter
plot_data_df <- plot_data_df[row_numbers_filter,]
## add a column by combining gene and phosphosite for filtering
plot_data_df <- plot_data_df %>%
  mutate(gene_phosphosite = paste0(Gene, "_", phosphosite))
## remove duplicates
plot_data_df <- plot_data_df[!duplicated(plot_data_df$gene_phosphosite),]
## transform data frame to matrix
plot_data_mat <- plot_data_df %>%
  select(tumor_bulk_aliquot_ids2plot, normal_bulk_aliquot_ids2plot)
rownames(plot_data_mat) <- plot_data_df$gene_phosphosite
## take out the reference intensity
plot_data_mat2 <- as.matrix(plot_data_mat) - as.vector(plot_data_df$ReferenceIntensity)

## make column annotation
## make column annotation
ca = HeatmapAnnotation(Sample_Type = ifelse(colnames(plot_data_mat2) %in% tumor_bulk_aliquot_ids2plot, "Tumor", "Normal"),
                       col = list(Sample_Type = c("Tumor" = "red", "Normal" = "green")))

## get color corresponding to values
col_protein = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

## plot heatmaps
p_mtor <- Heatmap(plot_data_mat2[(rowSums(!is.na(plot_data_mat2)) >= num_nonna),],
                  row_title = "mTORC1 Pathway",
                  name = "log2Intensity\n(Sample/Reference)", 
                  row_title_gp = gpar(fontsize = 12),
                  col = col_protein,
                  row_names_gp = gpar(fontsize = row_fontsize),
                  top_annotation = ca,
                  cluster_columns = F,
                  show_column_names = F,
                  cluster_rows = T)
p_mtor
## save heatmap to PNG
file2write <- paste0(dir_out, "Drug_Targets_Phosphoprotein.", run_id, ".png")
png(file2write, width = 2000, height = 400, res = 150)
print(p_mtor)
dev.off()


