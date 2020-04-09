# Yige Wu @ WashU 2020 Feb
## extrac RCC PDX from the merged gene expression table

# set up libraries and output directory -----------------------------------
## set up working directory and source functions and load libraries
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/ccRCC_Drug_shared.R")
## set run id 
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies --------------------------------------------
## input sample info to add in information about the model id, passage and TumorTissue info
sample_info_df <- readxl::read_excel(path = "./PDX-Pilot/DataFreeze/Sample_info/sampleInfo.v3/sampleInfo.washU_b1-b8.pdmr.other.passed.v3.20200130.from_jay.xlsx")
sample_info_df <- rbind(sample_info_batch8[,colnames(sample_info_df)], sample_info_df)
## input RNA expression
rna_exp_df <- fread("./Ding_Lab/Projects_Current/PDX-WashU/merged_GeneExp/washU.b1-b8.trans2gene.kallisto.tsv", data.table = F)
rna_exp_df %>% head()

# get the sample ids to filter-------------------------------------------
## get column names/analysis ids
analysis_ids <- colnames(rna_exp_df)
analysis_ids <- col_names[grepl(pattern = "resl", x = col_names, ignore.case = T)]
## only keep PDX, no PDXO or PDXS because the low sequencing depth from batch 8 QC
analysis_ids <- analysis_ids[analysis_ids %in% sample_info_df$Analysis_ID]
## get the sample ids
sample_ids <- plyr::mapvalues(x = analysis_ids, from = sample_info_df$Analysis_ID, to = as.vector(sample_info_df$SampleID))
sample_ids

# filter gene expression matrix -------------------------------------------
rna_exp_df_filtered <- rna_exp_df[,c("Gene", analysis_ids)]

# write output-------------------------------------------------------------------
## write output in analysis id as column names
write.table(x = rna_exp_df_filtered, file = paste0(dir_out, "RCC_PDX.GeneExp.TPM.AnalysisID.tsv"), quote = F, row.names = F, sep = "\t")
## write output in sample id as column names
rna_exp_df_filtered2 <- rna_exp_df_filtered
colnames(rna_exp_df_filtered2) <- c("Gene", sample_ids)
write.table(x = rna_exp_df_filtered2, file = paste0(dir_out, "RCC_PDX.GeneExp.TPM.SampleID.tsv"), quote = F, row.names = F, sep = "\t")


