# Yige Wu @ WashU 2022 Jan

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
files_process <- list.files(path = "./Resources/Analysis_Results/expression/protein/pathway/", recursive = T, pattern = "ORA_Results.tsv")
files_process
enricher_all_df <- NULL
for (file_tmp in files_process) {
  enricher_df <- fread(data.table = F, input = paste0("./Resources/Analysis_Results/expression/protein/pathway/", file_tmp))
  test_tmp <- str_split_fixed(string = file_tmp, pattern = "\\/", n = 3)[1]
  enricher_df$test <- test_tmp
  enricher_all_df <- rbind(enricher_df, enricher_all_df)
}
enricher_all_df <- enricher_all_df %>%
  mutate(number_genes_tested = str_split_fixed(string = GeneRatio, pattern = "\\/", n = 2)[,2]) %>%
  mutate(number_genes_tested = as.numeric(number_genes_tested)) %>%
  mutate(generatio_num = Count/number_genes_tested) %>%
  mutate(row_id = paste0(test, "|", ID))
# filter ------------------------------------------------------------------
# enricher_filtered_df <- enricher_all_df %>%
#   group_by(test) %>%
#   slice_min(pvalue, n = 5)

enricher_filtered_df <- enricher_all_df %>%
  filter(pvalue < 0.05)
enricher_filtered_df$Keep <- F
## step one: add one pathway at a time based on the max overlap
enricher_filtered_df$Keep[enricher_filtered_df$test == "ora_msigdb_H_CP_Cabo_human_protein_neg_markers_in_control" & enricher_filtered_df$ID %in% c("REACTOME_RNA_POLYMERASE_II_TRANSCRIPTION", "WP_LUNG_FIBROSIS", "PID_ERBB1_INTERNALIZATION_PATHWAY", "REACTOME_ION_CHANNEL_TRANSPORT", "REACTOME_PHASE_II_CONJUGATION_OF_COMPOUNDS")] <- T
enricher_filtered_df$Keep[enricher_filtered_df$test == "ora_msigdb_H_CP_Cabo_human_protein_pos_markers_in_control" & enricher_filtered_df$ID %in% c("REACTOME_ASPARAGINE_N_LINKED_GLYCOSYLATION", "REACTOME_SIGNALING_BY_WNT", "KEGG_LYSOSOME", "REACTOME_SARS_COV_INFECTIONS", "REACTOME_DISEASES_OF_METABOLISM")] <- T
enricher_filtered_df$Keep[enricher_filtered_df$test == "ora_msigdb_H_CP_Sap_human_protein_neg_markers_in_control" & enricher_filtered_df$ID %in% c("REACTOME_ADAPTIVE_IMMUNE_SYSTEM", "HALLMARK_G2M_CHECKPOINT", "REACTOME_KERATINIZATION", "WP_ERBB_SIGNALING_PATHWAY", "WP_IL18_SIGNALING_PATHWAY")] <- T
enricher_filtered_df$Keep[enricher_filtered_df$test == "ora_msigdb_H_CP_Sap_human_protein_pos_markers_in_control" & enricher_filtered_df$ID %in% c("NABA_ECM_AFFILIATED", "HALLMARK_XENOBIOTIC_METABOLISM", "HALLMARK_MITOTIC_SPINDLE", "REACTOME_EGR2_AND_SOX10_MEDIATED_INITIATION_OF_SCHWANN_CELL_MYELINATION", "REACTOME_MRNA_DECAY_BY_5_TO_3_EXORIBONUCLEASE")] <- T
enricher_filtered_df$Keep[enricher_filtered_df$test == "ora_msigdb_H_CP_CaboSap_human_protein_neg_markers_in_control" & enricher_filtered_df$ID %in% c("REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM", "REACTOME_ESR_MEDIATED_SIGNALING")] <- T
enricher_filtered_df$Keep[enricher_filtered_df$test == "ora_msigdb_H_CP_CaboSap_human_protein_pos_markers_in_control" & enricher_filtered_df$ID %in% c("REACTOME_CELLULAR_RESPONSES_TO_EXTERNAL_STIMULI", "REACTOME_HIV_INFECTION", "WP_MRNA_PROCESSING", "REACTOME_RAC1_GTPASE_CYCLE", "REACTOME_GOLGI_ASSOCIATED_VESICLE_BIOGENESIS")] <- T

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

## step 3: pick one pathway and go back to step 1 to add it
temp_df <- enricher_filtered_df %>%
  filter(test == "ora_msigdb_H_CP_CaboSap_human_protein_pos_markers_in_control") %>%
  # filter(test == "ora_msigdb_H_CP_CaboSap_human_protein_pos_markers_in_control") %>%
  arrange(desc(generatio_num))
enricher_top_df <- enricher_filtered_df %>%
  filter(Count >= 3) %>%
  filter(Keep)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ora_msigdb_H_CP_human_protein_markers_in_control.", run_id, ".tsv")
write.table(x = enricher_all_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "ora_msigdb_H_CP_human_protein_markers_in_control.top.", run_id, ".tsv")
write.table(x = enricher_top_df, file = file2write, quote = F, sep = "\t", row.names = F)
