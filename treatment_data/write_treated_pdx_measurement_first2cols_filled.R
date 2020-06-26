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

# input dependencies --------------------------------------------
## input treatment data
path_excel <- "./Resources/Treatment_Lists/RCC_PDX_Treatment_Lists/20200625_RCC_PDX_Treatment_MeasurementOnly.v1.xlsx"
tab_names <- excel_sheets(path = path_excel)
tab_names
list_all <- lapply(tab_names, function(x) read_excel(path = path_excel, sheet = x))
names(list_all) <- tab_names

# format ------------------------------------------------------------------
list_measurement_data_df <- list()
for (tab_tmp in tab_names) {
  ## extract current tab
  measurement_raw_df <- list_all[[tab_tmp]]
  
  ## divide into data and header
  measurement_header_df <- measurement_raw_df[1:3,]
  measurement_data_df <- measurement_raw_df[4:nrow(measurement_raw_df),]
  ## fill in header
  measurement_header_df <- as.data.frame(measurement_header_df)
  for (row_index in 1:2) {
    for (col_index in 4:ncol(measurement_header_df)) {
      value_now <- measurement_header_df[row_index, col_index]
      if (!is.na(value_now)) {
        next()
      } else {
        value_left <- measurement_header_df[row_index, col_index-1]
        measurement_header_df[row_index, col_index] <- value_left
      }
    }
  }
  ## make unique headers for the measurement data
  colnames_data_uniq <- sapply(1:ncol(measurement_header_df), function(i, data_df) { paste0(data_df[,i], collapse = "-") }, data_df = measurement_header_df)
  colnames_data_uniq[1:3] <- c("Treatment_Status", "Date", "Measurement_Type")
  colnames_data_uniq
  ## rename measurement data
  colnames(measurement_data_df) <- colnames_data_uniq
  measurement_data_df <- as.data.frame(measurement_data_df)
  ## fill in date and date annotation
  for (i in 2:nrow(measurement_data_df)) {
    value_now <- measurement_data_df$Date[i]
    if (!is.na(value_now)) {
      next()
    } else {
      value_previous <- measurement_data_df$Date[i-1]
      measurement_data_df$Date[i] <- value_previous
    }
  }
  for (i in 2:nrow(measurement_data_df)) {
    value_now <- measurement_data_df$Treatment_Status[i]
    if (!is.na(value_now)) {
      next()
    } else {
      value_previous <- measurement_data_df$Treatment_Status[i-1]
      measurement_data_df$Treatment_Status[i] <- value_previous
    }
  }
  file2write <- paste0(dir_out, "RCC_PDX.", tab_tmp, ".MeasurementOnly.", run_id, ".tsv")
  write.table(x = measurement_data_df, file = file2write, quote = F, sep = "\t", row.names = F)
  
  ## combine date annotation
  list_measurement_data_df[[tab_tmp]] <- measurement_data_df
}
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "RCC_PDX.", "MeasurementOnly.", run_id, ".xlsx")
write.xlsx(list_measurement_data_df, file = file2write)

