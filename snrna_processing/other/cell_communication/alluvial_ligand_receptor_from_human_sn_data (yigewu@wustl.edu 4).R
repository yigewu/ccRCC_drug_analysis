# Yige Wu @WashU Apr 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/plotting.R")
# install.packages("ggalluvial")
library(ggalluvial)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
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

# make plot data ----------------------------------------------------------
plotdata_df <- str_split_fixed(string = interacting_pair_plot_sorted, pattern = "_", n = 2) %>% as.data.frame()
colnames(plotdata_df) <- c("genesymbol_ligand", "genesymbol_receptor")
plotdata_df$genesymbol_ligand <- factor(x = plotdata_df$genesymbol_ligand, levels = unique(plotdata_df$genesymbol_ligand))
plotdata_df$genesymbol_receptor <- factor(x = plotdata_df$genesymbol_receptor, levels = unique(plotdata_df$genesymbol_receptor))

# plot --------------------------------------------------------------------
p <- ggplot(data = plotdata_df,
       aes(axis1 = genesymbol_ligand, axis2 = genesymbol_receptor)) +
  geom_alluvium(alpha = 0.7, color = "black") +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal()
p
file2write <- paste0(dir_out, "LigandGenes", ".png")
png(file2write, width = 800, height = 800, res = 150)
print(p)
dev.off()

