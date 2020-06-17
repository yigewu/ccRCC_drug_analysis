# Yige Wu @ WashU 2020 Jun

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
## input treatment data
path_excel <- "./Resources/Treatment_Lists/RCC_PDX_Treatment_Lists/20200615_RCC_PDX_Treatment_MeasurementOnly.v1.xlsx"
tab_names <- excel_sheets(path = path_excel)
tab_names
list_all <- lapply(tab_names, function(x) read_excel(path = path_excel, sheet = x))
names(list_all) <- tab_names

# format ------------------------------------------------------------------
measurement_long_all_df <- NULL
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
  ## fill in date
  for (i in 2:nrow(measurement_data_df)) {
    value_now <- measurement_data_df$Date[i]
    if (!is.na(value_now)) {
      next()
    } else {
      value_previous <- measurement_data_df$Date[i-1]
      measurement_data_df$Date[i] <- value_previous
    }
  }
  ## get length, width, volume in long data frame
  volume_wide_df <- measurement_data_df %>%
    filter(!is.na(Measurement_Type)) %>%
    filter(Measurement_Type %in% c("length/mm", "width/mm", "volume/mm3")) %>%
    select(-Treatment_Status)
  volume_long_df <- melt(data = volume_wide_df, id.vars = c("Date", "Measurement_Type"))
  volume_long_df <- volume_long_df %>%
    mutate(Group = str_split_fixed(string = variable, pattern = "-", n = 3)[,1]) %>%
    mutate(ID.Xiaolu = str_split_fixed(string = variable, pattern = "-", n = 3)[,2]) %>%
    mutate(TumorPiece = str_split_fixed(string = variable, pattern = "-", n = 3)[,3])
  colnames(volume_long_df)
  ## get body weight and dose in long data frame
  bw_dose_wide_df <- measurement_data_df %>%
    filter(!is.na(Measurement_Type)) %>%
    filter(Measurement_Type %in% c("B.W./g", "dose/uL", "PO/uL", "IP/uL", "other")) %>%
    select(-Treatment_Status)
  bw_dose_long_df <- melt(data = bw_dose_wide_df, id.vars = c("Date", "Measurement_Type"))
  bw_dose_long_df <- bw_dose_long_df %>%
    mutate(Group = str_split_fixed(string = variable, pattern = "-", n = 3)[,1]) %>%
    mutate(ID.Xiaolu = str_split_fixed(string = variable, pattern = "-", n = 3)[,2]) %>%
    mutate(TumorPiece = NA)
  
  ## get treatment status and comment for each date
  date_anno_df <- measurement_data_df %>%
    select(Date, Treatment_Status) %>%
    unique() %>%
    arrange(Treatment_Status)
  date_anno_df <- date_anno_df[!duplicated(date_anno_df$Date),]
  
  ## merge data
  measurement_long_df <- rbind(volume_long_df %>%
                                 filter(!is.na(value)) %>%
                                 select(Date, Group, ID.Xiaolu, TumorPiece, Measurement_Type, value),
                               bw_dose_long_df %>%
                                 filter(!is.na(value)) %>%
                                 select(Date, Group, ID.Xiaolu, TumorPiece, Measurement_Type, value))
  measurement_long_df <- merge(measurement_long_df, date_anno_df, by = c("Date"), all.x = T)
  measurement_long_df <- measurement_long_df %>%
    arrange(Date, Group, ID.Xiaolu, TumorPiece, Measurement_Type) %>%
    mutate(Tab = tab_tmp)
  ## add days from the start
  treatmenttime_start <- min(date_anno_df$Date[grepl(x = date_anno_df$Treatment_Status, pattern = "Treatment ON", ignore.case = T)])
  treatmenttime_start <- as.Date(treatmenttime_start)
  measurement_long_df <- measurement_long_df %>%
    mutate(Date = as.Date(Date)) %>%
    mutate(Date.TreatmentStart = treatmenttime_start) %>%
    mutate(Treatment.Days = Date - Date.TreatmentStart)
  ## add to the parent table
  measurement_long_all_df <- rbind(measurement_long_df, measurement_long_all_df)
}
unique(measurement_long_all_df$Treatment.Days)
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "RCC_PDX.", "Measurement.", "Long.", run_id, ".tsv")
write.table(x = measurement_long_all_df, file = file2write, sep = "\t", quote = F, row.names = F)

