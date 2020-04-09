# Yige Wu @ WashU 2018 Jan
## shared plotting functions for phospho_network

# library -----------------------------------------------------------------
packages = c(
  "ggplot2",
  "RColorBrewer",
  "pheatmap",
  "ggbeeswarm"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

# pheatmap ----------------------------------------------------------------

save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf <- function(x, filename, width=6, height=6) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}



