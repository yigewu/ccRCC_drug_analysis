packages = c(
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}