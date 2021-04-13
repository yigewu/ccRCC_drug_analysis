# Yige Wu @WashU Apr 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/plotting.R")
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the average expression
exp_wide_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snrna_processing/average_expression/unite_expression/unite_sct_data_humancells_mousecells_bycelltypeshorter/20210406.v1/SCT.data.AverageExpression.20210406.v1.tsv")

# identify genes to plot -------------------------------------------------
interacting_pair_plot_sorted <- c(
  ## tumor&stroma
  "VEGFA_FLT1", "VEGFA_KDR", "VEGFB_FLT1", "PGF_FLT1", "EFNA5_EPHA2", "EFNA5_EPHA4", 
  "EFNB2_EPHB1", "EFNB2_EPHA4", "EFNA1_EPHA7", "FGF2_FGFR2", "FGF2_FGFR1", "EFNA1_EPHA4", "NRG2_ERBB3", "FGF2_FGFR3",
  "HGF_MET", "FGF1_FGFR1",
  "EFNA5_EPHA3", "NRG3_ERBB4", "EFNA1_EPHA3", "NRG1_ERBB4", 
  "TGFA_EGFR", "EFNA5_EPHA7", "NRG1_ERBB3",
  ## tumor&tumor
  "TGFA_EGFR", "EFNA5_EPHA7", "EGF_EGFR", "EFNA5_EPHA4", "EFNA1_EPHA7", "NRG1_ERBB3", "FGF2_FGFR1", "NCAM1_FGFR1",
  ## stroma & stroma
  "VEGFA_FLT1", "PGF_FLT1", "EFNB2_EPHA4", "EFNB2_EPHB1", "EFNB2_EPHB4", "ANGPT2_TEK", "VEGFC_KDR", "VEGFC_FLT4", "EFNB1_EPHA4", "EFNA1_EPHA4", "FGF2_FGFR1", "EFNB1_EPHB4", "EFNA1_EPHA2",
  "VEGFB_FLT1", "ANGPT1_TEK", "FGF1_FGFR1",
  "EFNA5_EPHA7",
  "PDGFD_PDGFRB", "PDGFB_PDGFRB", "EFNA1_EPHA3", "EFNB2_EPHB6",
  "EFNA5_EPHA7", "TGFA_EGFR", 
  "EFNA5_EPHA3",
  "NRG3_ERBB4", "NTF3_NTRK3",
  ## tumor&immune
  "HBEGF_EGFR", "TGFA_EGFR", "IGF1_IGF1R", "HGF_MET",
  "IL34_CSF1R"
)
genes_plot <- str_split_fixed(string = interacting_pair_plot_sorted, pattern = "_", n = 2)[,1]
genes_plot <- unique(genes_plot)

# make plot data ----------------------------------------------------------
plotdata_wide_df <- exp_wide_df %>%
  filter(genesymbol_human %in% genes_plot)
plotdata_df <- melt(data = plotdata_wide_df)
summary(plotdata_df$value)
x_cap <- 4
plotdata_df <- plotdata_df %>%
  rename(cell_group = variable) %>%
  mutate(x_plot = ifelse(value > x_cap, x_cap, value))
plotdata_df$y_plot <- factor(x = plotdata_df$genesymbol_human, levels = rev(genes_plot))
## make colors
colors_cellgroup <- RColorBrewer::brewer.pal(n = 5, name = "Set1")
names(colors_cellgroup) <- unique(plotdata_df$cell_group)
# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot, color = cell_group), alpha = 0.7, size = 3)
p <- p + scale_color_manual(values = colors_cellgroup)
p <- p + theme_classic(base_size = 12)
p <- p + xlab("Normalized expression")
p <- p + theme(axis.text.y = element_text(size = 12), axis.title.y = element_blank())
p <- p + theme(axis.text.x = element_text(size = 12))
p <- p + scale_x_reverse() + scale_y_discrete(position = "right")
p <- p + theme(legend.position = "left")
file2write <- paste0(dir_out, "LigandGenes", ".png")
png(file2write, width = 800, height = 800, res = 150)
print(p)
dev.off()

