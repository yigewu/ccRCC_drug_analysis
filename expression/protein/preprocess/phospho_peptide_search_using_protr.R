# Yige Wu @ WashU 2021 Feb

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
packages_defaultinstall <- c("protr")
for (pkg_name_tmp in packages_defaultinstall) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  library(package = pkg_name_tmp, character.only = T)
}


# input dependencies ------------------------------------------------------


# get unique protein ids --------------------------------------------------
ids <- c("P00750", "P00751", "P00752")
getUniProt(id = ids)

