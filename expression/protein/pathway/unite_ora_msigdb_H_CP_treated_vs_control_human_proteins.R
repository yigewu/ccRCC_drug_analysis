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
files_process <- files_process[grepl(pattern = "20220209", x = files_process)]; files_process
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
### for each group, add the top pathway first
table(enricher_filtered_df$test)
enricher_filtered_df$Keep[enricher_filtered_df$test == "ora_msigdb_H_CP_Cabo_decreased_human_protein" &
                            enricher_filtered_df$ID %in% c("REACTOME_CELL_CYCLE", "REACTOME_TRANSLATION", "REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES",
                                                           "REACTOME_PROGRAMMED_CELL_DEATH", "HALLMARK_E2F_TARGETS", "HALLMARK_UV_RESPONSE_UP",
                                                           "REACTOME_CYTOPROTECTION_BY_HMOX1", "REACTOME_CELL_CELL_COMMUNICATION", "HALLMARK_PEROXISOME",
                                                           "REACTOME_VXPX_CARGO_TARGETING_TO_CILIUM")] <- T
# "REACTOME_DEVELOPMENTAL_BIOLOGY", "REACTOME_CELL_CYCLE", "REACTOME_TRANSLATION",
# "REACTOME_PROGRAMMED_CELL_DEATH", "HALLMARK_E2F_TARGETS", "HALLMARK_UV_RESPONSE_UP",
# "REACTOME_CYTOPROTECTION_BY_HMOX1", "REACTOME_CELL_CELL_COMMUNICATION", "HALLMARK_PEROXISOME",
# "REACTOME_VXPX_CARGO_TARGETING_TO_CILIUM", "BIOCARTA_ERAD_PATHWAY", "WP_MITOCHONDRIAL_CIV_ASSEMBLY",
# "KEGG_PORPHYRIN_AND_CHLOROPHYLL_METABOLISM"
enricher_filtered_df$Keep[enricher_filtered_df$test == "ora_msigdb_H_CP_Cabo_increased_human_protein" &
                            enricher_filtered_df$ID %in% c("NABA_MATRISOME", "HALLMARK_GLYCOLYSIS", "HALLMARK_APICAL_JUNCTION",
                                                           "HALLMARK_APOPTOSIS", "REACTOME_RAC1_GTPASE_CYCLE", "REACTOME_CILIUM_ASSEMBLY",
                                                           "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "REACTOME_L1CAM_INTERACTIONS", "REACTOME_TRANSCRIPTIONAL_ACTIVITY_OF_SMAD2_SMAD3_SMAD4_HETEROTRIMER",
                                                           "WP_ENVELOPE_PROTEINS_AND_THEIR_POTENTIAL_ROLES_IN_EDMD_PHYSIOPATHOLOGY")] <- T
enricher_filtered_df$Keep[enricher_filtered_df$test == "ora_msigdb_H_CP_Sap_decreased_human_protein" &
                            enricher_filtered_df$ID %in% c("REACTOME_TRANSLATION", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "REACTOME_TRANSPORT_OF_SMALL_MOLECULES", 
                                                           "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES", "REACTOME_PROTEIN_LOCALIZATION",
                                                           "HALLMARK_ADIPOGENESIS", "REACTOME_METABOLISM_OF_CARBOHYDRATES", "HALLMARK_GLYCOLYSIS",
                                                           "HALLMARK_E2F_TARGETS")] <- T
# , "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_MTORC1_SIGNALING",
# "REACTOME_FATTY_ACID_METABOLISM", "REACTOME_DISEASES_OF_METABOLISM", "REACTOME_COLLAGEN_FORMATION"
enricher_filtered_df$Keep[enricher_filtered_df$test == "ora_msigdb_H_CP_Sap_increased_human_protein" &
                            enricher_filtered_df$ID %in% c("REACTOME_VESICLE_MEDIATED_TRANSPORT", "REACTOME_SIGNALING_BY_RHO_GTPASES_MIRO_GTPASES_AND_RHOBTB3", "KEGG_FOCAL_ADHESION",
                                                           "NABA_MATRISOME", "KEGG_MAPK_SIGNALING_PATHWAY", "KEGG_ENDOCYTOSIS",
                                                           "REACTOME_SIGNALING_BY_GPCR", "HALLMARK_COMPLEMENT", "HALLMARK_HYPOXIA", 
                                                           "REACTOME_APOPTOSIS")] <- T
# "PID_PDGFRB_PATHWAY" , "HALLMARK_APOPTOSIS", "HALLMARK_MYOGENESIS",
# "REACTOME_SMOOTH_MUSCLE_CONTRACTION", "WP_MYOMETRIAL_RELAXATION_AND_CONTRACTION_PATHWAYS", "KEGG_DILATED_CARDIOMYOPATHY"
enricher_filtered_df$Keep[enricher_filtered_df$test == "ora_msigdb_H_CP_CaboSap_decreased_human_protein" &
                            enricher_filtered_df$ID %in% c("REACTOME_TRANSLATION", "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES", "HALLMARK_E2F_TARGETS",
                                                           "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "REACTOME_PROTEIN_LOCALIZATION", "HALLMARK_MYC_TARGETS_V1",
                                                           "HALLMARK_ADIPOGENESIS", "WP_CILIARY_LANDSCAPE", "REACTOME_CELL_CYCLE_CHECKPOINTS",
                                                           "HALLMARK_MTORC1_SIGNALING")] <- T
# , "REACTOME_FATTY_ACID_METABOLISM", "HALLMARK_FATTY_ACID_METABOLISM",
# "WP_NUCLEAR_RECEPTORS_METAPATHWAY", "REACTOME_CILIUM_ASSEMBLY", "HALLMARK_PROTEIN_SECRETION"
enricher_filtered_df$Keep[enricher_filtered_df$test == "ora_msigdb_H_CP_CaboSap_increased_human_protein" &
                            enricher_filtered_df$ID %in% c("NABA_MATRISOME", "WP_NUCLEAR_RECEPTORS_METAPATHWAY", "HALLMARK_GLYCOLYSIS",
                                                           "HALLMARK_HYPOXIA", "HALLMARK_APOPTOSIS", "KEGG_LYSOSOME", "REACTOME_CELLULAR_SENESCENCE",
                                                           "REACTOME_SYNTHESIS_OF_SUBSTRATES_IN_N_GLYCAN_BIOSYTHESIS", "REACTOME_MITOTIC_PROPHASE",
                                                           "REACTOME_TP53_REGULATES_METABOLIC_GENES")] <- T

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
  filter(test == "ora_msigdb_H_CP_Cabo_decreased_human_protein") %>%
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
pathwaynames_increased_process <- unique(enricher_filtered_df$ID[enricher_filtered_df$Keep & grepl(pattern = "increased", x = enricher_filtered_df$test)])
enricher_increasedprotein_df <- enricher_processed_df %>%
  filter(grepl(pattern = "increased", x = test)) %>%
  filter(ID %in% pathwaynames_increased_process)
pathwaynames_decreased_process <- unique(enricher_filtered_df$ID[enricher_filtered_df$Keep & grepl(pattern = "decreased", x = enricher_filtered_df$test)])
enricher_decreasedprotein_df <- enricher_processed_df %>%
  filter(grepl(pattern = "decreased", x = test)) %>%
  filter(ID %in% pathwaynames_decreased_process)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ora_msigdb_H_CP.treated_vs_control.human.proteins.", run_id, ".tsv")
write.table(x = enricher_processed_df, file = file2write, quote = F, sep = "\t", row.names = F)

file2write <- paste0(dir_out, "ora_msigdb_H_CP.treated_vs_control.human.proteins.increased.top.", run_id, ".tsv")
write.table(x = enricher_increasedprotein_df, file = file2write, quote = F, sep = "\t", row.names = F)

file2write <- paste0(dir_out, "ora_msigdb_H_CP.treated_vs_control.human.proteins.decreased.top.", run_id, ".tsv")
write.table(x = enricher_decreasedprotein_df, file = file2write, quote = F, sep = "\t", row.names = F)

