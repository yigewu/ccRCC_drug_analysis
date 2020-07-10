# Yige Wu @WashU Jun 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
source("./ccRCC_drug_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input KEGG
load("./Resources/Knowledge/Databases/2015-08-01_Gene_Set.RData")
## input HIF target genes
hiftarget_genes <- fread(data.table = F, input = "../ccRCC_snRNA/Resources/Analysis_Results/dependencies/write_hif_targets/20200428.v1/HIF_Target_Genes.20200428.v1.tsv")

# make data frame ---------------------------------------------------------
## add important pathways for RCC: https://www.genome.jp/kegg-bin/show_pathway?hsa05211+4233
### including HGF-MET pathway, PI3K, RAS-RAF-MEK-ERK
rccgenes_df1 <- data.frame(genesymbol = KEGG[["hsa05211\tRenal cell carcinoma"]], genesetname = "Renal cell carcinoma", source = "KEGG-hsa05211", category = "RCC Pathogenic Pathways")
rccgenes_df2 <- hiftarget_genes %>%
  rename(genesymbol = target_genesymbol) %>%
  mutate(genesetname = ifelse(source_genesymbol == "HIF1A", "HIF1A Targets", "EPAS1 Targets")) %>%
  mutate(source = "OmniPath") %>%
  mutate(category = "RCC Pathogenic Pathways") %>%
  select(genesymbol, genesetname, source, category)
rccgenes_df <- rbind(rccgenes_df1, rccgenes_df2)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "KEGG_RCC.", run_id, ".tsv")
write.table(x = rccgenes_df, file = file2write, quote = F, row.names = F, sep = "\t")

