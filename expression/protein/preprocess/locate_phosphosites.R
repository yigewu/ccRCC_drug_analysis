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
## input the phopho data
phospho_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Data_Files/Protein/01252021/LiDing_PDX_ccRCC_204phospho_library_Report.txt")
## input phosphopeptide location
peptide_loc_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/preprocess/phospho_peptide_search_using_predownloaded_fasta/20210203.v1/Phosphopeptide_Location.20210203.v1.tsv")

# preprocess --------------------------------------------------------------
phospho_filtered_df <- phospho_df %>%
  filter(grepl(pattern = "phospho", x = EG.ModifiedSequence, ignore.case = T))
phospho_filtered_uniq_df <- phospho_filtered_df %>%
  select(PG.ProteinGroups, PG.ProteinNames, PEP.StrippedSequence, EG.ModifiedSequence) %>%
  unique()
## split multiple accessions
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
proteinnames_vec <- sapply(phospho_filtered_uniq_df$PG.ProteinNames, function(x) {
  x_split <- str_split(string = x, pattern = "\\;")[[1]]
  return(x_split)
}, simplify = T)
proteinnames_vec <- unlist(proteinnames_vec)
phospho_filtered_uniq_acc_df <- data.frame(PG.ProteinGroups = rep(phospho_filtered_uniq_df$PG.ProteinGroups, idx_rep),
                                           PG.ProteinAccession = accessions_vec,
                                           PG.ProteinNames = rep(phospho_filtered_uniq_df$PG.ProteinNames, idx_rep),
                                           PG.ProteinName = proteinnames_vec,
                                           PEP.StrippedSequence = rep(phospho_filtered_uniq_df$PEP.StrippedSequence, idx_rep),
                                           EG.ModifiedSequence = rep(phospho_filtered_uniq_df$EG.ModifiedSequence, idx_rep))
## add peptide locatin
phospho_filtered_uniq_acc_df <- merge(x = phospho_filtered_uniq_acc_df, y = peptide_loc_df, 
                                      by = c("PG.ProteinGroups", "PG.ProteinAccession", "PEP.StrippedSequence"), all.x = T)

# map phophosite ----------------------------------------------------------
ptms_vec <- NULL
i <- 1
for (i in 1:nrow(phospho_filtered_uniq_acc_df)) {
  print(i)
  modified_peptide_tmp <- phospho_filtered_uniq_acc_df[i, "EG.ModifiedSequence"] %>% gsub(pattern = "_", replacement = "")
  unmodified_peptide_tmp <- phospho_filtered_uniq_acc_df[i, "PEP.StrippedSequence"] %>% as.vector()
  start_tmp <- phospho_filtered_uniq_acc_df[i, "start"] %>% as.vector()
  end_tmp <- phospho_filtered_uniq_acc_df[i, "end"] %>% as.vector()
  ## identify the positions and amino acids with motification string
  starts_mat <- str_locate_all(string = modified_peptide_tmp, pattern = "\\[")[[1]]
  ends_mat <- str_locate_all(string = modified_peptide_tmp, pattern = "\\]")[[1]]
  startends_df <- data.frame(start_offset = starts_mat[,"start"], end_offset = ends_mat[, "end"])
  startends_df <- startends_df %>%
    mutate(insertion_length = (end_offset - start_offset + 1)) %>%
    mutate(cumsum_insertion_length = cumsum(insertion_length)) %>%
    mutate(modified_offset = (end_offset - cumsum_insertion_length)) %>%
    mutate(modified_position = (start_tmp+modified_offset-1))
  startends_df$modified_aa <- sapply(startends_df$modified_offset, function(i, string) {
    chars_vec <- str_split(string = string, pattern = "")[[1]]
    if (i > 0) {
      char_i <- chars_vec[i]
    } else {
      char_i <- "Acetyl"
    }
    return(char_i)
  }, string = unmodified_peptide_tmp)
  startends_df <- startends_df %>%
    mutate(modification_name = paste0(modified_aa, modified_position))
  ## return string of the PTMs
  ptms_vec <- c(ptms_vec, paste0(startends_df$modification_name, collapse = "_"))
}
phospho_filtered_uniq_acc_df$PTM_Name <- ptms_vec

# ## identify all the modifications
# modified_peptide_split <- str_split(string = modified_peptide_tmp, pattern = "\\[|\\]")[[1]]
# ### extract elements at even indexes
# idxs_modification <- 2*(1:floor(length(modified_peptide_split)/2))
# modifications_tmp <- modified_peptide_split[idxs_modification]
# j <- 1
# modification_tmp <- modifications_tmp[j]; modification_noblack <- str_split(string = modification_tmp, pattern = " ")[[1]][1]
# 

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Phosphorylation_Sites.", run_id, ".tsv")
write.table(x = phospho_filtered_uniq_acc_df, file = file2write, quote = F, sep = "\t", row.names = F)

