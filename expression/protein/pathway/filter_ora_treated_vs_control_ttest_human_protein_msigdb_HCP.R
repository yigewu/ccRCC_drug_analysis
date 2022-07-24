# Yige Wu @ WashU 2022 March

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
  "data.table"
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


# input -------------------------------------------------------------------
enricher_all_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/pathway/ora_treated_vs_control_ttest_human_protein_msigdb_HCP/20220323.v1/Treated_vs_control.ORA_Results.tsv")

# process -----------------------------------------------------------------
enricher_all_df <- enricher_all_df %>%
  mutate(number_genes_tested = str_split_fixed(string = GeneRatio, pattern = "\\/", n = 2)[,2]) %>%
  mutate(number_genes_tested = as.numeric(number_genes_tested)) %>%
  mutate(generatio_num = Count/number_genes_tested)%>%
  mutate(test = paste0(group1, "_vs_control.", diff_direction)) %>%
  mutate(row_id = paste0(test, "|", ID))
enricher_filtered_df <- enricher_all_df %>%
  filter(pvalue < 0.05)
enricher_filtered_df$Keep <- F
## step one: add one pathway at a time based on the max overlap
### for each group, add the top pathway first and add one for each iteration
table(enricher_filtered_df$test)
enricher_filtered_df$Keep[enricher_filtered_df$test == "Cabo+ Sap_vs_control.down" &
                            enricher_filtered_df$ID %in% c("REACTOME_TRANSLATION", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "REACTOME_PROTEIN_LOCALIZATION",
                                                           "HALLMARK_E2F_TARGETS", "HALLMARK_BILE_ACID_METABOLISM")] <- T
enricher_filtered_df$Keep[enricher_filtered_df$test == "Cabo+ Sap_vs_control.up" &
                            enricher_filtered_df$ID %in% c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM", "HALLMARK_G2M_CHECKPOINT",
                                                           "REACTOME_NUCLEAR_ENVELOPE_BREAKDOWN", "WP_INTRACELLULAR_TRAFFICKING_PROTEINS_INVOLVED_IN_CMT_NEUROPATHY")] <- T
enricher_filtered_df$Keep[enricher_filtered_df$test == "Cabo_vs_control.down" &
                            enricher_filtered_df$ID %in% c("REACTOME_PROGRAMMED_CELL_DEATH", "REACTOME_M_PHASE", "REACTOME_FORMATION_OF_THE_CORNIFIED_ENVELOPE")] <- T
enricher_filtered_df$Keep[enricher_filtered_df$test == "Cabo_vs_control.up" &
                            enricher_filtered_df$ID %in% c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_INFLAMMATORY_RESPONSE", "REACTOME_RHOG_GTPASE_CYCLE",
                                                           "HALLMARK_IL2_STAT5_SIGNALING", "WP_UREA_CYCLE_AND_METABOLISM_OF_AMINO_GROUPS")] <- T
enricher_filtered_df$Keep[enricher_filtered_df$test == "Sap_vs_control.down" &
                            enricher_filtered_df$ID %in% c("REACTOME_ASPARAGINE_N_LINKED_GLYCOSYLATION", "REACTOME_TRANSPORT_OF_SMALL_MOLECULES", "REACTOME_CILIUM_ASSEMBLY",
                                                           "KEGG_OXIDATIVE_PHOSPHORYLATION", "REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES")] <- T
enricher_filtered_df$Keep[enricher_filtered_df$test == "Sap_vs_control.up" &
                            enricher_filtered_df$ID %in% c("HALLMARK_MITOTIC_SPINDLE", "HALLMARK_G2M_CHECKPOINT", "KEGG_FOCAL_ADHESION", 
                                                           "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_MYOGENESIS")] <- T

## step 2: run this
enricher_filtered_df$max_overlap_ratio <- sapply(enricher_filtered_df$row_id, function(row_id_tmp, test_df) {
  # sapply(head(enricher_filtered_df$row_id), function(row_id_tmp, test_df) {
  
  test_tmp <- test_df$test[test_df$row_id == row_id_tmp]
  generatio_tmp <- test_df$generatio_num[test_df$row_id == row_id_tmp]
  id_tmp <- test_df$ID[test_df$row_id == row_id_tmp]
  test_keep_df <- test_df[test_df$test == test_tmp & test_df$Keep & test_df$generatio_num >= generatio_tmp & test_df$ID != id_tmp,]
  
  genes_tmp <- str_split(string = test_df$geneID[test_df$row_id == row_id_tmp], pattern = "\\/")[[1]]
  max_overlap_ratio <- 0
  for (genestring_kept_tmp in test_keep_df$geneID) {
    genes_kept_tmp <- str_split(string = genestring_kept_tmp, pattern = "\\/")[[1]]
    genes_common_tmp <- intersect(genes_tmp, genes_kept_tmp)
    overlap_ratio_tmp <- length(genes_common_tmp)/length(genes_tmp)
    max_overlap_ratio <- max(c(max_overlap_ratio, overlap_ratio_tmp))
  }
  return(max_overlap_ratio)
}, test_df = enricher_filtered_df)

## step 3: pick one pathway and go back to step 1 to add one more pathway
temp_df <- enricher_filtered_df %>%
  filter(test == "Sap_vs_control.down") %>%
  arrange(desc(generatio_num))
View(temp_df)

# filter by all pathways chosen -------------------------------------------
enricher_top_df <- enricher_all_df[enricher_all_df$row_id %in% enricher_filtered_df$row_id[enricher_filtered_df$Keep],]
enricher_top_df <- enricher_top_df %>%
  group_by(test) %>%
  mutate(my_ranks = order(order(generatio_num, decreasing = T)))
enricher_nontop_df <- enricher_all_df[!(enricher_all_df$row_id %in% enricher_filtered_df$row_id[enricher_filtered_df$Keep]),]
enricher_nontop_df$my_ranks <- NA
enricher_processed_df <- rbind(enricher_top_df, enricher_nontop_df)
pathwaynames_increased_process <- unique(enricher_filtered_df$ID[enricher_filtered_df$Keep & enricher_filtered_df$diff_direction == "up"])
enricher_increasedprotein_df <- enricher_processed_df %>%
  filter(diff_direction == "up") %>%
  filter(ID %in% pathwaynames_increased_process)
pathwaynames_decreased_process <- unique(enricher_filtered_df$ID[enricher_filtered_df$Keep & enricher_filtered_df$diff_direction == "down"])
enricher_decreasedprotein_df <- enricher_processed_df %>%
  filter(diff_direction == "down") %>%
  filter(ID %in% pathwaynames_decreased_process)

# output ------------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_drug_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

file2write <- paste0(dir_out, "ora_msigdb_H_CP.treated_vs_control.human.proteins.", run_id, ".tsv")
write.table(x = enricher_processed_df, file = file2write, quote = F, sep = "\t", row.names = F)

file2write <- paste0(dir_out, "ora_msigdb_H_CP.treated_vs_control.human.proteins.increased.top.", run_id, ".tsv")
write.table(x = enricher_increasedprotein_df, file = file2write, quote = F, sep = "\t", row.names = F)

file2write <- paste0(dir_out, "ora_msigdb_H_CP.treated_vs_control.human.proteins.decreased.top.", run_id, ".tsv")
write.table(x = enricher_decreasedprotein_df, file = file2write, quote = F, sep = "\t", row.names = F)



