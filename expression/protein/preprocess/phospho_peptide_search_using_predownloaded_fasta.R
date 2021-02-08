# Yige Wu @ WashU 2021 Feb

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
packages_defaultinstall <- c("seqinr")
for (pkg_name_tmp in packages_defaultinstall) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the phopho data
phospho_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Data_Files/Protein/01252021/LiDing_PDX_ccRCC_204phospho_library_Report.txt")
## input FASTA data
fasta_list <- read.fasta(file = "./Resources/Knowledge/Databases/Uniprot/Phosphoprotein")

# preprocess --------------------------------------------------------------
phospho_filtered_df <- phospho_df %>%
  filter(grepl(pattern = "phospho", x = EG.ModifiedSequence, ignore.case = T))
phospho_filtered_uniq_df <- phospho_filtered_df %>%
  select(PG.ProteinGroups, PG.ProteinNames, PEP.StrippedSequence) %>%
  unique()
idx_rep <- sapply(phospho_filtered_uniq_df$PG.ProteinGroups, function(x) {
  x_split <- str_split(string = x, pattern = "\\;")[[1]]
  x_len <- length(x_split)
  return(x_len)
})
idx_rep <- as.numeric(idx_rep)
accessions_vec <- sapply(phospho_filtered_uniq_df$PG.ProteinGroups, function(x) {
  x_split <- str_split(string = x, pattern = "\\;")[[1]]
  return(x_split)
}, simplify = T)
accessions_vec <- unlist(accessions_vec)
## create unique accessions
phospho_filtered_uniq_acc_df <- data.frame(PG.ProteinGroups = rep(phospho_filtered_uniq_df$PG.ProteinGroups, idx_rep),
                                           PG.ProteinAccession = accessions_vec,
                                           PEP.StrippedSequence = rep(phospho_filtered_uniq_df$PEP.StrippedSequence, idx_rep))
## process fasta
fasta_names_df <- data.frame(Name = names(fasta_list))
fasta_names_df <- fasta_names_df %>%
  mutate(Accession_w_isoform = str_split_fixed(string = Name, pattern = "\\|", n = 3)[,2]) %>%
  mutate(Accession_Canonical = str_split_fixed(string = Accession_w_isoform, pattern = "\\-", n = 2)[,1]) %>%
  mutate(Protein_Name = str_split_fixed(string = Name, pattern = "\\|", n = 3)[,3])

# run by unique peptide ---------------------------------------------------
peptide_location_df <- NULL
for (i in 1:nrow(phospho_filtered_uniq_acc_df)) {
  print(i)
  peptide_tmp <- phospho_filtered_uniq_acc_df[i, "PEP.StrippedSequence"] %>% as.vector()
  accession_tmp <- phospho_filtered_uniq_acc_df[i, "PG.ProteinAccession"]
  canonical_fastaname_tmp <- fasta_names_df$Name[fasta_names_df$Accession_w_isoform == accession_tmp] %>% as.vector()
  protein_string_tmp <- getSequence(object = fasta_list[[canonical_fastaname_tmp]], as.string = T)[[1]]
  protein_string_tmp <- toupper(protein_string_tmp)
  peptide_location_tmp <- str_locate(string = protein_string_tmp, pattern = peptide_tmp)
  if (any(is.na(peptide_location_tmp))) {
    stop("meow")
  }
  peptide_location_df <- rbind(peptide_location_df, peptide_location_tmp)
}
View(peptide_location_df)
phospho_filtered_uniq_acc_df <- cbind(phospho_filtered_uniq_acc_df, peptide_location_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Phosphopeptide_Location.", run_id, ".tsv")
write.table(x = phospho_filtered_uniq_acc_df, file = file2write, quote = F, sep = "\t", row.names = F)

