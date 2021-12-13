# Yige Wu @ WashU 2021 Nov
## need to overlap with protein changes

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
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
## input the genes to plot
pathway2genes_df1 <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/rna/pathway/ora_msigdb_H_CP_Cabo_sensitive_markers/20211123.v1/ORA_Results.tsv")
pathway2genes_df2 <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/rna/pathway/ora_msigdb_H_CP_Cabo_resistant_markers/20211123.v1/ORA_Results.tsv")
pathway2genes_df3 <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/rna/pathway/ora_msigdb_H_CP_Sap_sensitive_markers/20211118.v1/ORA_Results.tsv")
pathway2genes_df4 <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/rna/pathway/ora_msigdb_H_CP_Sap_resistant_markers/20211118.v1/ORA_Results.tsv")
## input pathway pairwise comparison
# ora_obj1 <- readRDS(file = "./Resources/Analysis_Results/expression/rna/pathway/ora_msigdb_H_CP_Cabo_sensitive_markers/20211111.v1/ORA_Results.RDS")

# unite and select pathways -----------------------------------------------------------------
pathway2genes_df <- rbind(pathway2genes_df1 %>%
                            mutate(Comparison = "Cabo sensitive"),
                          pathway2genes_df2 %>%
                            mutate(Comparison = "Cabo resistant"),
                          pathway2genes_df3 %>%
                            mutate(Comparison = "Sap sensitive"),
                          pathway2genes_df4 %>%
                            mutate(Comparison = "Sap resistant"))
pathway2genes_sig_df <- pathway2genes_df[pathway2genes_df$p.adjust < 0.05,]
## select pathways
### select top 3 most significantly enriched pathways (jaccard similarity coefficient < 0.1) for each comparisons
pathway_selected <- c("REACTOME_TRANSLATION", "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
                      "REACTOME_HSF1_ACTIVATION", "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION", #"HALLMARK_WNT_BETA_CATENIN_SIGNALING",
                      "REACTOME_TRANSLATION", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", #"HALLMARK_FATTY_ACID_METABOLISM", 
                      "HALLMARK_MTORC1_SIGNALING",
                      "REACTOME_PCP_CE_PATHWAY", "HALLMARK_MYC_TARGETS_V1")

pathway2genes_filtered_df <- pathway2genes_df %>%
  filter(Description %in% pathway_selected) %>%
  group_by(Comparison) %>%
  mutate(rank_pvalue = rank(pvalue, ties.method = "first"))

# map each gene to pathway ------------------------------------------------
pathway2genes_list <- sapply(pathway2genes_filtered_df$geneID, function(x) {
  genes_vec <- str_split(string = x, pattern = "\\/")[[1]]
  return(genes_vec)
})
genes2filter <- unique(unlist(pathway2genes_list))
gene2pathway_df <- data.frame(GeneSet_Name = rep(x = pathway2genes_filtered_df$Description, sapply(X = pathway2genes_list, FUN = function(x) {
  length_vec <- length(x)
  return(length_vec)
})), GeneSymbol = unlist(pathway2genes_list))
gene2pathway_df <- unique(gene2pathway_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Cabo_Sap.Gene2TopPathway.", run_id, ".tsv")
write.table(x = gene2pathway_df, file = file2write, sep = "\t", row.names = F, quote = F)
file2write <- paste0(dir_out, "Cabo_Sap.Pathways2Genes.", run_id, ".tsv")
write.table(x = pathway2genes_df, file = file2write, sep = "\t", row.names = F, quote = F)
file2write <- paste0(dir_out, "Cabo_Sap.TopPathways2Genes.", run_id, ".tsv")
write.table(x = pathway2genes_filtered_df, file = file2write, sep = "\t", row.names = F, quote = F)
