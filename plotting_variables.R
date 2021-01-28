# make color palette for variant class ------------------------------------
cartocolors_df <- rcartocolor::cartocolors
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
