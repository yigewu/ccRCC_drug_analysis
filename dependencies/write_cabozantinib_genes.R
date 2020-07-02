# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# gather cabo genes --------------------------------------------------------
## Cabozantinib: VEGFR2 inhibitor, also inhibits c-Met, Ret, Kit, Flt-1/3/4, Tie2, and AXL
genes_cabo_targets <- c("MET", "AXL",
                        "KDR", "FLT3", "FLT1", "FLT4",
                        "KIT", 
                        "RET", ## https://www.ncbi.nlm.nih.gov/pubmed/23705946/
                        "NTRK2", "TEK",
                        "TIE2", "TYRO3", "TRKB",
                        "ROS1") ## https://pubmed.ncbi.nlm.nih.gov/25351743/
genes_hgf_met <- c("HGF", "MET", 
                   "GAS6", "AXL", ## https://pubmed.ncbi.nlm.nih.gov/32055714/
                   "SRC") ## https://www.futuremedicine.com/doi/full/10.2217/fon-2019-0021?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub++0pubmed&
## Growth arrest-specific 6/AXL forms a complex with the proto-oncogene tyrosine kinase inhibitor (SRC) and activates the MET receptor in an HGF-independent manner
genes_vegf <- c("FLT1", "KDR", "FLT3", "FLT4", 
                "NRP1", "NRP2",
                "VEGFA", "VEGFB", "VEGFC", "VEGFD", "VEGFE")
cabo_substrate_genes <- c("CYP3A4", "MRP2")
genes_cabo_relevant <- c(genes_cabo_targets,
                         genes_hgf_met, genes_vegf, 
                         cabo_substrate_genes)
genes_cabo_relevant <- unique(genes_cabo_relevant)

s# make data frame ---------------------------------------------------------
cabogenes_df <- data.frame(genesymbol = genes_cabo_relevant, pathwaname = "Cabozantinib", pathwayid = NA, pathwaysource = "Literature Review")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Cabozantinib_Genes.", run_id, ".tsv")
write.table(x = genes_df, file = file2write, quote = F, row.names = F, sep = "\t")
