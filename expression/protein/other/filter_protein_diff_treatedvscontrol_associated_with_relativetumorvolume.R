# Yige Wu @ WashU 2022 Feb

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
# ## set run id
# version_tmp <- 1
# run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
# ## set output directory
# dir_out <- paste0(makeOutDir(), run_id, "/")
# dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input pathway enrichment results
ora_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/pathway/unite_ora_msigdb_H_CP_treated_vs_control_human_proteins/20220216.v1/ora_msigdb_H_CP.treated_vs_control.human.proteins.20220216.v1.tsv")
## input correlation result
cor_result_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/test/cor_1month_proteinschanges_vs_relativetumorvolumevscontrol_bytreatment/20220217.v1/spearman.1month.controlproteins_vs_relativetumorvolume.20220217.v1.tsv")
## input the differnetial protein test result
exp_diff_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/test/wilcox_paired_diff_protein_treated_vs_control/20220209.v1/Wilcox.Paired.Each_Treated_Group_vs_CT.20220209.v1.tsv")

# filter for combo --------------------------------------------------------
unique(ora_df$test)
unique(cor_result_df$treatment)
test_tmp <- "ora_msigdb_H_CP_CaboSap_decreased_human_protein"
test_tmp <- "ora_msigdb_H_CP_CaboSap_increased_human_protein"

path_id_tmp <- "REACTOME_TRANSLATION"
path_id_tmp <- "HALLMARK_MTORC1_SIGNALING"
path_id_tmp <- "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES"
path_id_tmp <- "HALLMARK_MYC_TARGETS_V1"
path_id_tmp <- "KEGG_LYSOSOME"
path_id_tmp <- "REACTOME_CELLULAR_SENESCENCE"
path_id_tmp <- "NABA_MATRISOME"
path_id_tmp <- "WP_NUCLEAR_RECEPTORS_METAPATHWAY"

treatment_tmp <- "Cabo+ Sap"
genestring_tmp <- ora_df$geneID[ora_df$test == test_tmp & ora_df$ID == path_id_tmp]
genesymbols_tmp <- str_split(string = genestring_tmp, pattern = "\\/")[[1]]
cor_tmp_df <- cor_result_df %>%
  mutate(genesymbol = str_split_fixed(string = ID, pattern = "_", n = 2)[,1]) %>%
  filter(treatment == treatment_tmp) %>%
  filter(genesymbol %in% genesymbols_tmp) %>%
  arrange(p.value)

genestring_tmp <- ora_df$geneID[ora_df$test == test_tmp & !is.na(ora_df$my_ranks)]
genesymbols_tmp <- unlist(str_split(string = genestring_tmp, pattern = "\\/"))
cor_tmp_df <- cor_result_df %>%
  mutate(genesymbol = str_split_fixed(string = ID, pattern = "_", n = 2)[,1]) %>%
  filter(treatment == treatment_tmp) %>%
  filter(genesymbol %in% genesymbols_tmp) %>%
  arrange(p.value)

ora_tmp_df <- ora_df %>%
  filter(test == test_tmp)


exp_diff_tmp_df <- exp_diff_df %>%
  filter(group1 == treatment_tmp) %>%
  filter(PG.Genes %in% genesymbols_tmp) %>%
  # filter(PG.Genes %in% exp_diff_df$PG.Genes[exp_diff_df$pvalue < 0.05 & exp_diff_df$diff_estimate < 0 & exp_diff_df$group1 == "Sap"]) %>%
  arrange(diff_estimate)

# examine Cabo ------------------------------------------------------------
test_tmp <- "ora_msigdb_H_CP_Cabo_increased_human_protein"
path_id_tmp <- "NABA_MATRISOME"
treatment_tmp <- "Cabo"

genestring_tmp <- ora_df$geneID[ora_df$test == test_tmp & !is.na(ora_df$my_ranks)]
genesymbols_tmp <- unlist(str_split(string = genestring_tmp, pattern = "\\/"))
cor_tmp_df <- cor_result_df %>%
  mutate(genesymbol = str_split_fixed(string = ID, pattern = "_", n = 2)[,1]) %>%
  filter(treatment == treatment_tmp) %>%
  filter(genesymbol %in% genesymbols_tmp) %>%
  arrange(p.value)


genestring_tmp <- ora_df$geneID[ora_df$test == test_tmp & ora_df$ID == path_id_tmp]
genesymbols_tmp <- str_split(string = genestring_tmp, pattern = "\\/")[[1]]
cor_tmp_df <- cor_result_df %>%
  mutate(genesymbol = str_split_fixed(string = ID, pattern = "_", n = 2)[,1]) %>%
  filter(treatment == treatment_tmp) %>%
  filter(genesymbol %in% genesymbols_tmp) %>%
  arrange(p.value)
exp_diff_tmp_df <- exp_diff_df %>%
  filter(group1 == treatment_tmp) %>%
  filter(PG.Genes %in% genesymbols_tmp) %>%
  arrange(diff_estimate) %>%
  unique()


# examine Cabo ------------------------------------------------------------
test_tmp <- "ora_msigdb_H_CP_Sap_increased_human_protein"
path_id_tmp <- "NABA_MATRISOME"
treatment_tmp <- "Sap"
genestring_tmp <- ora_df$geneID[ora_df$test == test_tmp & ora_df$ID == path_id_tmp]
genesymbols_tmp <- str_split(string = genestring_tmp, pattern = "\\/")[[1]]
cor_tmp_df <- cor_result_df %>%
  mutate(genesymbol = str_split_fixed(string = ID, pattern = "_", n = 2)[,1]) %>%
  filter(treatment == treatment_tmp) %>%
  filter(genesymbol %in% genesymbols_tmp) %>%
  arrange(p.value)
exp_diff_tmp_df <- exp_diff_df %>%
  filter(group1 == treatment_tmp) %>%
  filter(PG.Genes %in% genesymbols_tmp) %>%
  arrange(diff_estimate) %>%
  unique()

genestring_tmp <- ora_df$geneID[ora_df$test == test_tmp & !is.na(ora_df$my_ranks)]
genesymbols_tmp <- unlist(str_split(string = genestring_tmp, pattern = "\\/"))
cor_tmp_df <- cor_result_df %>%
  mutate(genesymbol = str_split_fixed(string = ID, pattern = "_", n = 2)[,1]) %>%
  filter(treatment == treatment_tmp) %>%
  filter(genesymbol %in% genesymbols_tmp) %>%
  arrange(p.value)


