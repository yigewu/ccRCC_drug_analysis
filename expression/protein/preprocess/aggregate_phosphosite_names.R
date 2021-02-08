# Yige Wu @ WashU 2021 Feb

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
## input phosphosite location
ptm_names_byprotein_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/preprocess/locate_phosphosites/20210204.v1/Phosphorylation_Sites.20210204.v1.tsv")

# aggregate by human and mice ---------------------------------------------
ptm_names_byprotein_df <- ptm_names_byprotein_df %>%
  mutate(Protein_Species = ifelse(grepl(pattern = "Human", x = PG.ProteinName, ignore.case = T), "Human", "Mouse"))
ptm_names_byproteingroup_df <- ptm_names_byprotein_df %>%
  group_by(PG.ProteinGroups, PG.ProteinNames, EG.ModifiedSequence) %>%
  summarise(PTM_Name_Aggr = paste0(unique(PTM_Name), collapse = "|"),
            PTM_Name_Human = paste0(unique(PTM_Name[Protein_Species == "Human"]), collapse = "|"), 
            PTM_Name_Mouse = paste0(unique(PTM_Name[Protein_Species == "Mouse"]), collapse = "|"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Phosphorylation_Sites.Aggregated.", run_id, ".tsv")
write.table(x = ptm_names_byproteingroup_df, file = file2write, quote = F, sep = "\t", row.names = F)

