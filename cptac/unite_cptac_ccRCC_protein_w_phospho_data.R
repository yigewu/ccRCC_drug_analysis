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
  "data.table"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
source("./ccRCC_drug_analysis//functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input data --------------------------------------------------------------
## input phospho data
pho_df <- fread("../../CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/phosphoproteome/6_CPTAC3_CCRCC_Phospho_abundance_phosphosite_protNorm=2_CB.tsv", data.table = F)
## input protein data
protein_df <- fread("../../CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/proteome/6_CPTAC3_CCRCC_Whole_abundance_gene_protNorm=2_CB.tsv", data.table = F)
## input bulk meta data
bulk_meta_tab <- fread("../../CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/cptac-metadata.csv")

# process aliquot ids --------------------------------------------
### get the ids for the normal samples
normal_bulk_aliquot_ids2plot <- bulk_meta_tab$Specimen.Label[bulk_meta_tab$Set.A == "yes" & bulk_meta_tab$Type == "Normal"]
normal_bulk_aliquot_ids2plot
case_ids2plot <- mapvalues(x = normal_bulk_aliquot_ids2plot, from = bulk_meta_tab$Specimen.Label, to = bulk_meta_tab$Case.ID)
case_ids2plot
tumor_bulk_aliquot_ids2plot <- mapvalues(x = case_ids2plot, from = bulk_meta_tab$Case.ID[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"], to = as.vector(bulk_meta_tab$Specimen.Label[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"]))
tumor_bulk_aliquot_ids2plot


# preprocess phospho data -----------------------------------------------------------------
header_col_names_keep <- c("Protein_id", "Gene", "Phosphosite", "Index", "ReferenceIntensity")
pho_df2 <- pho_df %>%
  mutate(Phosphosite = str_split_fixed(string = Index, pattern = "_", n = 7)[,7]) %>%
  filter(Phosphosite != "")  %>%
  mutate(Protein_id = paste0(Gene, "_", Phosphosite)) %>%
  select(header_col_names_keep, tumor_bulk_aliquot_ids2plot, normal_bulk_aliquot_ids2plot)

# preprocess protein data -------------------------------------------------
protein_df2 <- protein_df %>%
  mutate(Protein_id = paste0(Index, "_Protein")) %>%
  mutate(Phosphosite = "") %>%
  rename(Gene = Index) %>%
  rename(Index = Proteins) %>%
  select(header_col_names_keep, tumor_bulk_aliquot_ids2plot, normal_bulk_aliquot_ids2plot)

# combine -----------------------------------------------------------------
combined_df <- rbind(pho_df2, protein_df2)
## make normalized data
combined_data_norm_mat <- as.matrix(combined_df[,c(tumor_bulk_aliquot_ids2plot, normal_bulk_aliquot_ids2plot)]) - as.vector(combined_df$ReferenceIntensity)
combined_data_norm_df <- cbind(combined_df[,header_col_names_keep], combined_data_norm_mat)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Protein_Phosphorylation_combined.", run_id, ".tsv")
write.table(x = combined_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "Protein_Phosphorylation_combined.refnormalized.", run_id, ".tsv")
write.table(x = combined_data_norm_df, file = file2write, quote = F, sep = "\t", row.names = F)
