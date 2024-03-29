# set up libraries and output directory -----------------------------------
## getting the path to the current script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
path_this_script <- thisFile()
## set working directory
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input -------------------------------------------------------------------
dirs_input <- c("/diskmnt/primary/ccRCC_Drug/")
md5_vec <- NULL
filename_vec <- NULL
path_vec <- NULL
for (dir_input in dirs_input) {
  files_input <- list.files(path = dir_input, recursive = T)
  files_input <- files_input[grepl(pattern = "md5", x = files_input)]
  # print(files_input)
  for (file_tmp in files_input) {
    print(file_tmp)
    text_tmp <- fread(data.table = F, input = paste0(dir_input, file_tmp))
    md5_vec[file_tmp] <- colnames(text_tmp)[1]
    filename_tmp <- str_split(string = file_tmp, pattern = "\\/")[[1]]
    filename_tmp = filename_tmp[length(filename_tmp)]
    filename_vec[file_tmp] <- filename_tmp
    path_vec[file_tmp] <- paste0(dir_input, file_tmp)
  }
}
md5_df <- data.frame(md5sum = md5_vec, file_name = filename_vec, file_path = path_vec)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC.PDX.snRNA.md5sum.", run_id, ".tsv")
write.table(x = md5_df, file = file2write, sep = "\t", row.names = F, quote = F)
