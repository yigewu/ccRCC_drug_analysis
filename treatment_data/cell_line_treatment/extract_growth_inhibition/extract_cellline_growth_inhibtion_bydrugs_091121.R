# Yige Wu @ WashU 2021 Aug

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
library(readxl)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies --------------------------------------------
dir_excel <- "./ccRCC_Alamar_Blue/ccRCC20210911-rep3/rep 3- 20210911-calculations/"

## input drug and concentration mapping
drug_conc_map_df <- read_xlsx(path = "./ccRCC_Alamar_Blue/ccRCC20210815-rep1/ccRCC drug conc plan 08092021.xlsx", range = "A93:G98", sheet = "experimental plan rep 2", col_names = T)

# process -----------------------------------------------------------------
filenames_excel <- list.files(path = dir_excel, recursive = T)
filenames_excel <- filenames_excel[!grepl(x = filenames_excel, pattern = "\\~")]
gi_group_sup_df <- NULL

filename_tmp <- filenames_excel[1]
for (filename_tmp in filenames_excel) {
  path_tmp <- paste0(dir_excel, filename_tmp)
  if (grepl(x = filename_tmp, pattern = "293T and RCC4 rep 2")) {
    cellline_name <- str_split_fixed(string = filename_tmp, pattern = "_", n = 6)[5]
    if(cellline_name == "RCC") {
      cellline_name <- "RCC-4-VHL-"
    }
  } else {
    cellline_name <- str_split(string = filename_tmp, pattern = " ", n = 2)[[1]][1]
  }
  
  gi_df <- read_xlsx(path = path_tmp, range = "B81:M86", sheet = "raw data", col_names = F)
  gi_vec <- as.vector(t(gi_df))
  group_df <- read_xlsx(path = path_tmp, range = "B89:M94", sheet = "raw data", col_names = F)
  group_vec <- as.vector(t(group_df))
  gi_group_df <- data.frame(growth_inhibition = gi_vec/100,
                            groupid = group_vec)
  gi_group_df <- gi_group_df %>%
    filter(!is.na(groupid)) %>%
    group_by(groupid) %>%
    mutate(repid = row_number()) %>%
    mutate(cellline = cellline_name) %>%
    mutate(drug_token = str_split_fixed(string = groupid, pattern = "\\/", n = 2)[,1]) %>%
    mutate(drug_token = gsub(pattern = " ", replacement = "", x = drug_token)) %>%
    mutate(conc_token = str_split_fixed(string = groupid, pattern = "\\/", n = 2)[,2]) %>%
    mutate(conc_token = gsub(pattern = " ", replacement = "", x = conc_token))
  
  gi_group_sup_df <- rbind(gi_group_sup_df, gi_group_df)
}
table(gi_group_sup_df$cellline)
# map drug and concentration ----------------------------------------------
drug_conc_map_df$drug_token <- paste0("drug", rownames(drug_conc_map_df))
drug_conc_long_df <- melt(data = drug_conc_map_df, id.vars = c("Drug", "Targt Gene", "unit", "drug_token"))
drug_conc_long_df <- drug_conc_long_df %>%
  mutate(conc_token = gsub(pattern = " ", replacement = "", x = variable)) %>%
  rename(conc_value = value)
gi_group_annotated_df <- merge(x = gi_group_sup_df, 
                         y = drug_conc_long_df %>%
                           select(Drug, drug_token, conc_token, conc_value), by = c("drug_token", "conc_token"), all.x = T)
gi_group_annotated_df <- gi_group_annotated_df %>%
  filter()
which(is.na(gi_group_annotated_df$cellline))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "RCC_Cellline.", "Growth_Inhibition.", run_id, ".tsv")
write.table(x = gi_group_annotated_df, file = file2write, sep = "\t", quote = F, row.names = F)

