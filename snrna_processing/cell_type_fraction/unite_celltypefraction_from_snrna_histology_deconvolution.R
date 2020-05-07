# Yige Wu @WashU Map 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input snRNA based cell type fraction
celltypefrac_snrna_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/cell_type_fraction/calculate_celltypefraction_from_4sampleintegration/20200501.v1/CellTypeFraction_from_4sampleintegration.20200501.v1.tsv")
## input histology-based cell type fraction
celltypefrac_hist_df <- readxl::read_xlsx(path = "./Resources/H&E_Images/HE_Review_CellTypeFraction.20200430.v1.xlsx")
# ## input deconvolution
# celltypescore_bulkrna_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/xCell/xCell_RCC_PDX.GeneExp.TPM.SampleID_xCell_0836042420.txt")
# ## input metadata
# idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/write_sequencing_status/20200316.v1/RCC_PDX_Related_Samples.Sequencing_Status.20200316.v1.tsv")

# format snRNA-based ------------------------------------------------------
## select tumor and vasculature from snRNA-based cell type fraction
celltypefrac_snrna_selected_df <- celltypefrac_snrna_df %>%
  filter(Cell_Type.Short %in% c("Tumor cells", "Endothelial cells"))
celltypefrac_snrna_selected_wide_df <- dcast(data = celltypefrac_snrna_selected_df, formula = Id_Model + Id_Sample ~ Cell_Type.Short, value.var = "Fraction_CellType_Sample")
celltypefrac_snrna_selected_wide_df <- celltypefrac_snrna_selected_wide_df %>%
  rename(Frac.Endothelial.snRNA = `Endothelial cells`) %>%
  rename(Frac.Tumorcells.snRNA = `Tumor cells`)

# merge -------------------------------------------------------------------
## merge histology based cell type fraction with snRNA-based
celltypefrac_merged_df <- merge(celltypefrac_snrna_selected_wide_df, celltypefrac_hist_df, by = c("Id_Sample"), all.x = T)
## arrange
celltypefrac_merged_df <- celltypefrac_merged_df %>%
  select(Id_Model, Id_Sample, Frac.Tumorcells.snRNA, Frac.Endothelial.snRNA, Frac.Tumorcells.Histology, Frac.Endothelial.Histology)
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "CellTypeFraction_from_snRNA_Histology.", run_id, ".tsv")
write.table(x = celltypefrac_merged_df, file = file2write, sep = "\t", quote = F, row.names = F)
