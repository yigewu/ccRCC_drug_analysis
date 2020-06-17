# Yige Wu @ WashU 2020 Feb
## merge PDX mutation + CNV status

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

## input mutation matrix
mut_mat <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/Resources/Analysis_Results/mutation/generate_bulk_mutation_table/20200215.v1/RCC_PDX.Mutation_Matrix.20200215.v1.tsv", data.table = F)

## input CNV matrix
cnv_mat <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/Resources/Analysis_Results/copy_number/make_arm_level_copy_number_matrix/20200215.v1/RCC_PDX.Arm_Level_CNV_Matrix.20200215.v1.tsv", data.table = F) 

# merge data --------------------------------------------------------------
## prepare mutation data for merging
mut_mat_merge <- melt(data = mut_mat, id.vars = "Hugo_Symbol")
mut_mat_merge <- mut_mat_merge %>%
  dplyr::rename(Info = Hugo_Symbol)
mut_mat_merge$SampleID <- sapply(mut_mat_merge$variable, function(v) {
  sampleid <- sample_info_df$SampleID[grepl(pattern = v, x = sample_info_df$Analysis_ID)]
  return(sampleid)
})

## prepare CNV data for merging
cnv_mat_merge <- melt(data = cnv_mat, id.vars = "Chromosome_Arm")
cnv_mat_merge <- cnv_mat_merge %>%
  dplyr::rename(Info = Chromosome_Arm)
cnv_mat_merge$SampleID <- sapply(cnv_mat_merge$variable, function(v) {
  sampleid <- sample_info_df$SampleID[grepl(pattern = v, x = sample_info_df$Analysis_ID)]
  return(sampleid)
})

## merge mutation and CNVs
merge_df <- rbind(mut_mat_merge, cnv_mat_merge)
merge_mat1 <- dcast(data = merge_df, formula = Info ~ SampleID, value.var = "value")
### sort
rownames(merge_mat1) <- merge_mat1$Info
merge_mat1 <- merge_mat1[c("VHL", "SETD2", "BAP1", "PBRM1", "KDM5C", "TP53", "3p", "5q", "14q"),]
merge_mat1[is.na(merge_mat1)] <- ""

## prepare ModelID and NCI_Passage info ready to merge
col_names <- data.frame(SampleID = c("SampleID", colnames(merge_mat1)[-1]),
                        ModelID = c("ModelID", plyr::mapvalues(x = colnames(merge_mat1)[-1], from = sample_info_df$SampleID, to = sample_info_df$ModelID)),
                        NCI_Passage = c("NCI_Passage", plyr::mapvalues(x = colnames(merge_mat1)[-1], from = sample_info_df$SampleID, to = sample_info_df$NCI_Passage)),
                        TumorTissue = c("TumorTissue", plyr::mapvalues(x = colnames(merge_mat1)[-1], from = sample_info_df$SampleID, to = sample_info_df$TumorTissue)))
col_names_t <- data.frame(t(col_names))
colnames(col_names_t) <- colnames(merge_mat1)

## merge all
merge_mat <- rbind(col_names_t, merge_mat1)

## sort 
colnames_sort <- col_names %>%
  mutate(ModelID.No = as.numeric(gsub(x = ModelID, pattern = "RESL", replacement = ""))) %>%
  arrange(ModelID.No, NCI_Passage)
colnames_sort <- c("Info", as.vector(colnames_sort$SampleID)[-nrow(colnames_sort)])
colnames_sort
merge_mat <- merge_mat[,colnames_sort]

# write output ------------------------------------------------------------
write.table(x = merge_mat, file = paste0(dir_out, "RCC_PDX.Somatic_Variant_Merged.", run_id, ".tsv"), quote = F, row.names = F, sep = "\t")

