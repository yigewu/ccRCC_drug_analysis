# Yige Wu @ WashU 2018 Jan
## shared plotting functions for phospho_network

# library -----------------------------------------------------------------
packages = c(
  "ggplot2",
  "RColorBrewer",
  "rcartocolor",
  "ggrepel"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
cartocolors_df <- rcartocolor::cartocolors

# make color palette for variant class ------------------------------------
# rcartocolor::display_carto_all()
cartocolors_temps <- cartocolors_df[cartocolors_df$Name == "Temps", "n7"][[1]]
cartocolors_tropic <- cartocolors_df[cartocolors_df$Name == "Tropic", "n7"][[1]]
variant_class_colors <- c(cartocolors_temps[1:4], 
                          cartocolors_tropic[4], 
                          cartocolors_temps[c(5,6,7)],
                          "white")
names(variant_class_colors) <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", 'Splice_Site', 
                                 "Silent", 
                                 "Missense_Mutation", "In_Frame_Ins", "In_Frame_Del",
                                 "None")


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



