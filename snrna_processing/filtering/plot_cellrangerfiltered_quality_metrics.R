# Yige Wu @WashU Apr 2020
## plot quality metrics for cell ranger filtered cells

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/load_pkgs.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/functions.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/ccRCC_drug_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies --------------------------------------------
## input the paths to the quality metrics
path_qm_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/Resources/Analysis_Results/quality_control/filtering/make_quality_metrics/20200401.v1/path_quality_metrics_per_barcode.20200401.v1.tsv", data.table = F)

# process each sample and make plots --------------------------------------
for (sample_id in path_qm_df$sample_id) {
  path_qm <- path_qm_df$path_output[path_qm_df$sample_id == sample_id]
  qm_df <- fread(input = path_qm, data.table = F)
  
  ## make scatterplot plot for number of UMI counts per species per barcode
  plot_data_df <- qm_df %>%
    filter(!is.na(barcode_raw)) %>%
    select(barcode_raw, `GRCh38-3.0.0.premrna`, mm10_premrna, call) %>%
    mutate(x = `GRCh38-3.0.0.premrna`) %>%
    mutate(y = mm10_premrna) %>%
    mutate(barcode_species = ifelse(call == "GRCh38-3.0.0.premrna", "Human",
                                    ifelse(call == "mm10_premrna", "Mouse", "Multiplet")))
  
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df, mapping = aes(x = x, y = y, color = barcode_species), alpha = 0.5)
  p <- p + coord_fixed(ratio = 1)
  p <- p + xlab(label = "nCount_RNA Mapped to\nHuman_Reference")
  p <- p + ylab(label = "nCount_RNA Mapped to\nMouse_Reference")
  p <- p + ggtitle(label = paste0(sample_id, ": ", nrow(plot_data_df), " Cell-Ranger-Filtered Barcodes"),
                   subtitle = paste0(length(which(plot_data_df$call == "GRCh38-3.0.0.premrna")), " Human Cells; ",
                                     length(which(plot_data_df$call == "mm10_premrna")), " Mouse Cells; ",
                                     length(which(plot_data_df$call == "Multiplet")), " Multiplet Cells; "))
  p <- p + theme(legend.position = "bottom")
  file2write <- paste0(dir_out, "scatterplot_ncount_rna_human_vs_mouse.", sample_id, ".", run_id, ".png")
  png(filename = file2write, width = 1500, height = 500, res = 150)
  print(p)
  dev.off()
  
  ## make violin plot for mitoRatio per species per barcode
  plot_data_df <- qm_df %>%
    filter(!is.na(barcode_raw)) %>%
    select(barcode_raw, mitoRatio, call) %>%
    mutate(barcode_species = ifelse(call == "GRCh38-3.0.0.premrna", "Human",
                                    ifelse(call == "mm10_premrna", "Mouse", "Multiplet"))) %>%
    mutate(x = barcode_species) %>%
    mutate(y = mitoRatio)
  p <- ggplot(data = plot_data_df, mapping = aes(x = x, y = y, fill = barcode_species, color = barcode_species))
  p <- p + geom_violin(alpha = 0.5)
  p <- p + geom_point(alpha = 0.2)
  p <- p + geom_quasirandom(alpha = 0.2, width = 0.2)
  p <- p + xlab(label = "Barcode Species")
  p <- p + ylab(label = "Ratio of Mitochondrial Gene Read Counts Per Barcode")
  p <- p + ggtitle(label = paste0(sample_id, ": ", nrow(plot_data_df), " Cell-Ranger-Filtered Barcodes"),
                   subtitle = paste0(length(which(plot_data_df$call == "GRCh38-3.0.0.premrna")), " Human Cells; ",
                                     length(which(plot_data_df$call == "mm10_premrna")), " Mouse Cells; ",
                                     length(which(plot_data_df$call == "Multiplet")), " Multiplet Cells; "))
  p <- p + theme(title = element_text(size = 8))
  p
  file2write <- paste0(dir_out, "violin_mitoRatio.", sample_id, ".", run_id, ".png")
  png(filename = file2write, width = 600, height = 800, res = 150)
  print(p)
  dev.off()
  
  ## make violin plot for nCount_RNA per species per barcode
  plot_data_df <- qm_df %>%
    filter(!is.na(barcode_raw)) %>%
    select(barcode_raw, nCount_RNA, call) %>%
    mutate(barcode_species = ifelse(call == "GRCh38-3.0.0.premrna", "Human",
                                    ifelse(call == "mm10_premrna", "Mouse", "Multiplet"))) %>%
    mutate(x = barcode_species) %>%
    mutate(y = nCount_RNA)
  p <- ggplot(data = plot_data_df, mapping = aes(x = x, y = y, fill = barcode_species))
  p <- p + geom_violin(alpha = 0.5)
  p <- p + geom_point(alpha = 0.2)
  p <- p + geom_quasirandom(alpha = 0.2, width = 0.2)
  p <- p + xlab(label = "Barcode Species")
  p <- p + ylab(label = "Number of RNA Reads  Per Barcode")
  p <- p + ggtitle(label = paste0(sample_id, ": ", nrow(plot_data_df), " Cell-Ranger-Filtered Barcodes"),
                   subtitle = paste0(length(which(plot_data_df$call == "GRCh38-3.0.0.premrna")), " Human Cells; ",
                                     length(which(plot_data_df$call == "mm10_premrna")), " Mouse Cells; ",
                                     length(which(plot_data_df$call == "Multiplet")), " Multiplet Cells; "))
  p <- p + theme(title = element_text(size = 8))
  p
  file2write <- paste0(dir_out, "violin_nCount_RNA.", sample_id, ".", run_id, ".png")
  png(filename = file2write, width = 600, height = 800, res = 150)
  print(p)
  dev.off()
  
  ## make violin plot for nFeature_RNA per species per barcode
  plot_data_df <- qm_df %>%
    filter(!is.na(barcode_raw)) %>%
    select(barcode_raw, nFeature_RNA, call) %>%
    mutate(barcode_species = ifelse(call == "GRCh38-3.0.0.premrna", "Human",
                                    ifelse(call == "mm10_premrna", "Mouse", "Multiplet"))) %>%
    mutate(x = barcode_species) %>%
    mutate(y = nFeature_RNA)
  p <- ggplot(data = plot_data_df, mapping = aes(x = x, y = y, fill = barcode_species))
  p <- p + geom_violin(alpha = 0.5)
  p <- p + geom_point()
  p <- p + geom_quasirandom(alpha = 0.1, width = 0.2)
  p <- p + xlab(label = "Barcode Species")
  p <- p + ylab(label = "Number of Genes Per Barcode")
  p <- p + ggtitle(label = paste0(sample_id, ": ", nrow(plot_data_df), " Cell-Ranger-Filtered Barcodes"),
                   subtitle = paste0(length(which(plot_data_df$call == "GRCh38-3.0.0.premrna")), " Human Cells; ",
                                     length(which(plot_data_df$call == "mm10_premrna")), " Mouse Cells; ",
                                     length(which(plot_data_df$call == "Multiplet")), " Multiplet Cells; "))
  p <- p + theme(title = element_text(size = 8))
  p
  file2write <- paste0(dir_out, "violin_nFeature_RNA.", sample_id, ".", run_id, ".png")
  png(filename = file2write, width = 600, height = 800, res = 150)
  print(p)
  dev.off()
}

