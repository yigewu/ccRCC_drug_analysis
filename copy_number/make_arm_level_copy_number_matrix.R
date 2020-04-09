# Yige Wu @ WashU 2020 Feb
## annotate PDX CNV status

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
## input the gistic arm-level data
gistic2_arm_cnv_df1 <- fread(input = "./PDX-Pilot/DataFreeze/CopyNumber/v2.gistic2.seg5kb.v20200114/tumorOnly.broad_values_by_arm.txt", data.table = F)
gistic2_arm_cnv_df2 <- fread(input = "./PDX-Pilot/DataFreeze/CopyNumber/v2.gistic2.seg5kb.v20200114/tumorNor.broad_values_by_arm.txt", data.table = F)

# filter results ----------------------------------------------------------
## set chromosome arms to filter
chr_arms2keep <- c("3p", "5q", "14q")
## filter by RESL and chromosome arms
gistic2_arm_cnv_df1_filtered <- gistic2_arm_cnv_df1[gistic2_arm_cnv_df1$`Chromosome Arm` %in% chr_arms2keep, grepl(pattern = "RESL", x = colnames(gistic2_arm_cnv_df1))]
gistic2_arm_cnv_df2_filtered <- gistic2_arm_cnv_df2[gistic2_arm_cnv_df2$`Chromosome Arm` %in% chr_arms2keep, grepl(pattern = "RESL", x = colnames(gistic2_arm_cnv_df2))]
## merge data
gistic2_arm_cnv_df <- data.frame(Chromosome_Arm = gistic2_arm_cnv_df1$`Chromosome Arm`[gistic2_arm_cnv_df1$`Chromosome Arm` %in% chr_arms2keep])
gistic2_arm_cnv_df <- cbind(gistic2_arm_cnv_df, cbind(gistic2_arm_cnv_df1_filtered, gistic2_arm_cnv_df2_filtered))

# write output ------------------------------------------------------------
write.table(x = gistic2_arm_cnv_df, file = paste0(dir_out, "RCC_PDX.Arm_Level_CNV_Matrix.", run_id, ".tsv"), quote = F, row.names = F, sep = "\t")


