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
## input gene-to-pathway
gene2pathway_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/rna/pathway/unite_ora_msigdb_H_CP_Cabo_Sap_markers/20211117.v1/Cabo_Sap.Gene2TopPathway.20211117.v1.tsv")
## input pbrm1 daps
degs_cabo_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/rna/get_foldchange/filter_Cabo_sensitive_resistant_genes_overlapping_proteins/20211110.v2/Cabo_related_genes.20211110.v2.tsv")
degs_sap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/rna/get_foldchange/filter_Sap_sensitive_resistant_genes_overlapping_proteins/20211110.v2/Sap_related_genes.20211110.v2.tsv")
## input pathway enrichment results
ora_df <- fread(data.table = F, input = "./Resources/Analysis_Results/expression/rna/pathway/unite_ora_msigdb_H_CP_Cabo_Sap_markers/20211117.v1/Cabo_Sap.TopPathways2Genes.20211117.v1.tsv")

# make plot data ----------------------------------------------------------
colnames_merge <- intersect(colnames(degs_cabo_df), colnames(degs_sap_df))
plotdata_source_df <- rbind(degs_cabo_df[,colnames_merge] %>%
                              mutate(DAP_group = paste0("PBRM1_", DAP_direction)),
                            degs_sap_df[,colnames_merge] %>%
                              mutate(DAP_group = paste0("BAP1_", DAP_direction)))
plotdata_source_df <- merge(x = plotdata_source_df, y = gene2pathway_df, by.x = c("Gene"), by.y = c("GeneSymbol"), all.x = T)  
## 
plotdata_df <- plotdata_source_df %>%
  filter(!is.na(GeneSet_Name)) %>%
  group_by(DAP_group) %>%
  select(GeneSet_Name) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 0)
plotdata_df <- merge(x = plotdata_df, 
                     y = ora_df %>%
                       select(Description, p.adjust, Comparison, rank_pvalue),
                     by.x = c("GeneSet_Name", "DAP_group"), by.y = c("Description", "Comparison"), all.x = T)
plotdata_df <- plotdata_df %>%
  mutate(log10FDR = -log10(p.adjust))
## order the pathways
pathway_selected <- c("HALLMARK_MTORC1_SIGNALING", "WP_EGFEGFR_SIGNALING_PATHWAY","REACTOME_RHO_GTPASE_CYCLE", "WP_DNA_REPAIR_PATHWAYS_FULL_NETWORK",
                      "WP_FOCAL_ADHESION", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "PID_EPHA_FWDPATHWAY", "WP_TGFB_SIGNALING_IN_THYROID_CELLS_FOR_EPITHELIALMESENCHYMAL_TRANSITION",
                      "WP_FOCAL_ADHESION", "REACTOME_RAC1_GTPASE_CYCLE", "HALLMARK_HYPOXIA")
pathway_ordered <- unique(pathway_selected)
plotdata_df$GeneSet_Name <- factor(x = plotdata_df$GeneSet_Name, levels = rev(pathway_ordered))
## make color patette
color_pal <- RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")
## renmake pathway labels
pathway_label_df <- data.frame(label = c("mTORC1 signaling", "EGF-EGFR signaling", "Rho GTPase cycle", "DNA repair", 
                                         "Focal adhesion", "TGFA signaling via NFKB", "EPHA forward signaling", "TGFB signaling for EMT", 
                                         "Rac1 GTPase cycle", "Hypoxia"),
                               pathway_name = c("HALLMARK_MTORC1_SIGNALING", "WP_EGFEGFR_SIGNALING_PATHWAY","REACTOME_RHO_GTPASE_CYCLE", "WP_DNA_REPAIR_PATHWAYS_FULL_NETWORK",
                                                "WP_FOCAL_ADHESION", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "PID_EPHA_FWDPATHWAY", "WP_TGFB_SIGNALING_IN_THYROID_CELLS_FOR_EPITHELIALMESENCHYMAL_TRANSITION",
                                                "REACTOME_RAC1_GTPASE_CYCLE", "HALLMARK_HYPOXIA"))

# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = DAP_group, y = GeneSet_Name, size = Freq, color = rank_pvalue), shape = 16)
p <- p + scale_size_area()
p <- p + scale_color_gradientn(colours = rev(color_pal), breaks = c(1, 3, 5, 7, 9))
p <- p + theme_classic()
p <- p + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, color = "black"),
               axis.title = element_blank(),
               axis.text.y = element_text(color = "black"),
               legend.position = "right", legend.direction = "horizontal")
p <- p + guides(color = guide_legend(title = "Rank by\np-value", nrow = 2, title.position = "top"),
                size = guide_legend(title = "Number of\npeaks", nrow = 2, title.position = "top"))
p <- p + scale_y_discrete(breaks = pathway_label_df$pathway_name, labels = pathway_label_df$label)
file2write <- paste0(dir_out, "bubble.colorByRankpvalue.legend.pdf")
pdf(file2write, width = 3.85, height = 2.75, useDingbats = F)
print(p)
dev.off()


# plot by log10FDR ----------------------------------------------------------------
# p <- ggplot()
# p <- p + geom_point(data = plotdata_df, mapping = aes(x = DAP_group, y = GeneSet_Name, size = Freq))
# p <- p + scale_size_area()
# file2write <- paste0(dir_out, "bubble.plain.png")
# png(file2write, width = 1000, height = 800, res = 150)
# print(p)
# dev.off()
# 
# p <- ggplot()
# p <- p + geom_point(data = plotdata_df, mapping = aes(x = DAP_group, y = GeneSet_Name, size = Freq, fill = log10FDR), color = "black", shape = 21)
# p <- p + scale_size_area()
# # p <- p + scale_color_gradient(low = "white smoke", high = "red")
# p <- p + scale_fill_gradientn(colours = color_pal)
# p <- p + theme_classic()
# p
# file2write <- paste0(dir_out, "bubble.colorBylog10FDR.pdf")
# pdf(file2write, width = 6, height = 4, useDingbats = F)
# print(p)
# dev.off()
# file2write <- paste0(dir_out, "bubble.colorlog10FDR.png")
# png(file2write, width = 1000, height = 800, res = 150)
# print(p)
# dev.off()

