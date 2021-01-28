# Yige Wu @ WashU 2020 Jun

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
library(readxl)
library(openxlsx)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
dir_files <- "./Resources/Treatment_Lists/RCC_PDX_Treatment_Lists/Treatment_List_By_Cohort/"
filenames <- list.files(path = dir_files)
filenames

# extract volumes-----------------------------------------------------------------
filename_tmp <- filenames[2]
## input measurements
measure_df <- read_excel(path = paste0(dir_files, filename_tmp), sheet = "MeasurementsAll", na = "NA")
volumes_df <- measure_df %>%
  filter(Measurement_Type == "volume/mm3")

# calculate relative volumes ----------------------------------------------
## remove before treatment data
idxs_treatmenton <- which(volumes_df$Treatment_Status == "Treatment on")
## transform data frame values to numeric
volumes_num_mat <- sapply(volumes_df[min(idxs_treatmenton):nrow(volumes_df),4:ncol(volumes_df)], as.numeric)
## replace 0s
volumes_num_mat[volumes_num_mat == 0] <- NA
## remove column values after the first NA value
volumes_num_clean_mat <- apply(volumes_num_mat, 2, function(x) {
  idx_nas <- which(is.na(x))
  x_new <- x
  
  if (length(idx_nas) > 0) {
    x_new[min(idx_nas):length(x_new)] <- NA
  }
  return(x_new)
})

# remove unwanted data ----------------------------------------------------
if (grepl(x = filename_tmp, pattern = "RESL10_B1")) {
  volumes_num_clean_mat[, "CT-12477-right"] <- NA
  volumes_num_clean_mat[, "CT-12481-right"] <- NA
  volumes_num_clean_mat[, "CT-12481-left"] <- NA
  volumes_num_clean_mat[, "Cab-12479-right"] <- NA
  volumes_num_clean_mat[, "Sap-12472-right"] <- NA
  volumes_num_clean_mat[, "Cab+Sap-12464-right"] <- NA
  volumes_num_clean_mat[, "Cab+Sap-12474-left"] <- NA
}
volumes_clean_df <- cbind(volumes_df[min(idxs_treatmenton):nrow(volumes_df),1:3], volumes_num_clean_mat)

## calculate relative tumor volumn compared to the begining of the treatment
relative_volume_mat <- sweep(x = volumes_num_clean_mat, 2, as.vector(volumes_num_clean_mat[1,]), `/`)*100
relative_volume_df <- cbind(volumes_df[min(idxs_treatmenton):nrow(volumes_df),1:3], relative_volume_mat)

# calculate average relative volumes --------------------------------------
relative_volume_long_df <- reshape2::melt(data = relative_volume_df, id.vars = c("Treatment_Status", "Date", "Measurement_Type"))
relative_volume_long_df <- relative_volume_long_df %>%
  mutate(Treatment_Group = str_split_fixed(string = variable, pattern = "-", n = 3)[,1])
avg_relative_volume_long_df <- relative_volume_long_df %>%
  dplyr::group_by(Treatment_Status, Measurement_Type, Treatment_Group, Date) %>%
  summarise(Avg_Relative_Volume = mean(x = value, na.rm = T))
avg_relative_volume_df <- dcast(data = avg_relative_volume_long_df, formula = Treatment_Status+Date+Measurement_Type ~ Treatment_Group, value.var = "Avg_Relative_Volume")
avg_relative_volume_df <- avg_relative_volume_df %>%
  arrange(Date)

# calculate relative volume standard deviation ----------------------------
sd_relative_volume_long_df <- relative_volume_long_df %>%
  dplyr::group_by(Treatment_Status, Measurement_Type, Treatment_Group, Date) %>%
  summarise(SD_Relative_Volume = sd(x = value, na.rm = T))
sd_relative_volume_df <- dcast(data = sd_relative_volume_long_df, formula = Treatment_Status+Date+Measurement_Type ~ Treatment_Group, value.var = "SD_Relative_Volume")
sd_relative_volume_df <- sd_relative_volume_df %>%
  arrange(Date)

# write output ------------------------------------------------------------
## make output directory
data_name <- gsub(x = filename_tmp, pattern = ".xlsx", replacement = "")
dir_out_now <- paste0(dir_out, data_name, "/")
dir.create(path = dir_out_now)
## write outputs
file2write <- paste0(dir_out_now, data_name, ".",  "TumorVolume", ".", run_id, ".tsv")
write.table(x = volumes_clean_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out_now, data_name, ".",  "RelativeTumorVolume", ".", run_id, ".tsv")
write.table(x = relative_volume_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out_now, data_name, ".",  "RelativeTumoreVolume.Average.ByTreatmentGroup", ".", run_id, ".tsv")
write.table(x = avg_relative_volume_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out_now, data_name, ".",  "RelativeTumorVolume.SD.ByTreatmentGroup", ".", run_id, ".tsv")
write.table(x = sd_relative_volume_df, file = file2write, quote = F, sep = "\t", row.names = F)

