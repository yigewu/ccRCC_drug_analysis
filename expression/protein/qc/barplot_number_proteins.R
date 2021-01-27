# Yige Wu @ WashU 2021 Jan

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
setwd(dir_base)
source("./ccRCC_drug_analysis/load_pkgs.R")
source("./ccRCC_drug_analysis/functions.R")
source("./ccRCC_drug_analysis/variables.R")
source("./ccRCC_drug_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the protein data
protein_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/protein/preprocess/normalize_DIA_protein_data/20210111.v1/RCC_PDX.DIA_Protein.Log2.QuantileNormalized.20210111.v1.tsv")

# count protein number ----------------------------------------------------
plot_data_df <- data.frame(sample_id = colnames(protein_df[,5:ncol(protein_df)]),
                           number_proteins = colSums(!is.na(protein_df[,5:ncol(protein_df)])))
median(x = plot_data_df$number_proteins)

# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_col(data = plot_data_df, mapping = aes(x = sample_id, y = number_proteins))
p <- p + theme_classic()
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
               axis.text.y = element_text(size = 12))
p

# write  ------------------------------------------------------------------
file2write <- paste0(dir_out, "Number_Proteins.png")
png(file2write, width = 2000, height = 800, res = 150)
print(p)
dev.off()
