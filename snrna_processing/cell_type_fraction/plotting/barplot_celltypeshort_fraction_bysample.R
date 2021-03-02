# Yige Wu @WashU Mar 2021
## plot dimplot with cell type annotated

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input fraction
celltype_frac_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/cell_type_fraction/calculate_celltypefraction_bysample/20210225.v1/CellTypeFraction.20210225.v1.tsv")

# prepare plot data -------------------------------------------------------
plot_data_df <-celltype_frac_df %>%
  mutate(Cell_Group_Plot = Cell_Type.Short) %>%
  mutate(Treatment_Group = str_split_fixed(string = Id_Sample, pattern = "-", n = 3)[,3]) %>%
  mutate(Treatment_Group = gsub(x = Treatment_Group, pattern = "2", replacement = "")) %>%
  mutate(Id_Model = str_split_fixed(string = Id_Sample, pattern = "-", n = 3)[,1])
plot_data_df$Treatment_Group <- factor(x = plot_data_df$Treatment_Group, levels = c("CT", "Cabo", "Sap", "Cabo_Sap"))
plot_data_df$Id_Sample <- factor(plot_data_df$Id_Sample)

# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_col(data = plot_data_df, mapping = aes(x = Treatment_Group, y = Fraction_CellType_Sample, fill = Cell_Group_Plot), position = "stack")
p <- p + facet_grid(Id_Model~. , scales = "free", space = "free", drop = T)
# p <- p + scale_x_discrete(breaks = plot_data_df$Id_Sample, labels = plot_data_df$Treatment_Group)
p <- p + coord_flip()
# p <- p + scale_fill_manual(values = colors_cellgroup13)
file2write <- paste0(dir_out, "Cell_Group_Plot_Composition.", "png")
png(file2write, width = 800, height = 400, res = 150)
print(p)
dev.off()

