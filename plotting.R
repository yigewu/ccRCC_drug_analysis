# Yige Wu @ WashU 2018 Jan
## shared plotting functions for phospho_network

# library -----------------------------------------------------------------
packages = c(
  "ggplot2",
  "RColorBrewer",
  "ggrepel"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
