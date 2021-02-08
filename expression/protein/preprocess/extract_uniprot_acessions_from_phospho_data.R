# Yige Wu @ WashU 2021 Jan

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
## input the phopho data
phospho_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Data_Files/Protein/01252021/LiDing_PDX_ccRCC_204phospho_library_Report.txt")

# get unique accessions ---------------------------------------------------
uniprotid_string_uniq <- unique(phospho_df$PG.ProteinAccessions)
length(uniprotid_string_uniq)
uniprotid_uniq <- sapply(uniprotid_string_uniq, function(x) { str_split(string = x, pattern = ";")[[1]]})
uniprotid_uniq <- unique(unlist(uniprotid_uniq))
length(uniprotid_uniq)
uniprotid_uniq_df <- data.frame(Acession = uniprotid_uniq)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Uniprot_Accessions.", run_id, ".tsv")
write.table(x = uniprotid_uniq_df, file = file2write, sep = "\t", quote = F, row.names = F)
