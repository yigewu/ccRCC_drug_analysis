# Yige Wu @ WashU 2020 May

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

# input dependencies --------------------------------------------
## input all ccRCC samples data info
data_status_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_status/write_bulk_data_status/20200612.v1/RCC_PDX_Related_Samples.PDXNet_B1_9.HTAN_B1.Data_Status.20200612.v1.tsv")
## input mutations of SMGs in wide data frame
# mut_wide_df <- fread(data.table = F, input = "./Resources/Analysis_Results/mutation/extract_smg_mutations/20200526.v1/RCC_PDX.somaticMut.SMG.Wide.20200526.v1.tsv")
mut_wide_df <- fread(data.table = F, input = "./Resources/Analysis_Results/mutation/extract_smg_mutations/20200528.v1/RCC_PDX.somaticMut.SMG.Wide.20200528.v1.tsv")
## input mutations of SMGs in wide data frame
cnv_wide_df <- fread(data.table = F, input = "./Resources/Analysis_Results/copy_number/extract_cnvgenes_gene_level_cnv/20200526.v1/RCC_PDX.Gene_Level_CNV.Known_CNV_Genes.Wide.CN.20200526.v1.tsv")

# change column names -----------------------------------------------------
## change mutation column names
colnames(mut_wide_df)
colnames(mut_wide_df) <- c("Analysis_ID", paste0("Mut.", colnames(mut_wide_df)[-1]))
## change CNV column names
colnames(cnv_wide_df)
colnames(cnv_wide_df) <- c("Analysis_ID", paste0("CN.", colnames(cnv_wide_df)[-1]))

# merge -------------------------------------------------------------------
mut_cnv_df <- merge(data_status_df %>%
                      filter(!is.na(Analysis_ID.WES)) %>%
                      filter(Group != "Human_Normal"), 
                    mut_wide_df, by.x = c("Analysis_ID.WES"), by.y = c("Analysis_ID"), all.x = T)
mut_cnv_df <- merge(mut_cnv_df, 
                    cnv_wide_df, by.x = c("Analysis_ID.WES"), by.y = c("Analysis_ID"), all.x = T)

# fill in NAs -------------------------------------------------------------
colnames_colanno <- colnames(mut_cnv_df)[grepl(pattern = "Mut", x = colnames(mut_cnv_df))]
index_na <- is.na(mut_cnv_df[,colnames_colanno])
mut_cnv_df[,colnames_colanno][index_na] <- "None"

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "RCC_PDX.", "SMG_Mutation.", "Known_Genes_CNV.", run_id, ".tsv")
write.table(x = mut_cnv_df, file = file2write, quote = F, sep = "\t", row.names = F)
