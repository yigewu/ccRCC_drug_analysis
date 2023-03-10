## set fixed parameters for Seurat
num_pcs <- 40
num_var_features <- 3000
findclusters_res <- 0.5

## set prefix for feature names for human genes and mouse genes
prefix_human_feature <- "GRCh38-3.0.0.premrna-"
prefix_mouse_feature <- "mm10-premrna---------"

# gene lists --------------------------------------------------------------
## SMG source 1: Sato et al, Nature Genetics, 2013
## SMG source 1: TCGA, Nature, 2013
ccRCC_SMGs <- c("VHL", "PBRM1", "BAP1", "TCEB1", "SETD2", "TP53", "PIK3CA", "MTOR", "KDM5C", "PTEN","TSC1",
                "CCNB2",
                "MSR1", "TXNIP", "BTNL3", "SLITRK6", "RHEB", "ARID1A", "NPNT", 
                "FPGT", "MUDENG", "TET2", "MUC4", "MLLT10", "MSGN1", "KRT32", "M6PR", "RPL14", "GRB7", "CSMD3", "DNHD1", "NLRP12", "VMO1", "OR4C13", "KCNMA1", "LMAN2L", "ZNF536", "YIPF3")
ccRCC_drivers <- c("VHL", "PBRM1", "BAP1", "TCEB1", "SETD2", "TP53", "PIK3CA", "MTOR", "KDM5C", "PTEN","TSC1")


# Two examples
# of cullin–RING ubiquitin ligase system molecular assemblies using CUL2 or CUL5 (left) and CUL3 (right) that interact with the BC-box protein–
# Elongin C–Elongin B complex and BTB protein, respectively, to recruit substrate for ubiquitination and subsequent degradation. VHL and KEAP1 are examples of BC-box and BTB proteins, respectively, that recruit HIF and NRF2 proteins for ubiquitin-mediated degradation. 
vhl_complex_genes <- c("VHL", 
                       "TCEB1", "TCEB2",
                       "ELOC", "ELOB",
                       "RBX1", "RNF7",
                       "CUL2", "CUL5")
keap1_nrf2_complex_genes <- c("NFE2L2", "KEAP1", "CUL3", "RBX1")
ubiquitin_proteasome_genes <- c(vhl_complex_genes,
                                keap1_nrf2_complex_genes,
                                "USP24", "NEDD4", "WWP2", "URB1", "USP34")
swisnf_genes <- c("PBRM1", "ARID1A", "ARID1B", "ARID2", "SMARCA2", "SMARCA4", "SMARCB1", "SMARCC1", "SMARCC2", "SMARCD2", "SMARCD2")
other_epigeneticregulator_genes <- c("SETD2", 
                                     "BAP1", 
                                     "KDM5C", "KDM6A",
                                     "TET2")
pi3k_mtor_genes <- c("EGFR", "ERBB3", "FGFR3", "FGFR4", "IGF1R", 
                     "PIK3CA", "PIK3CB", "PIK3CG","PTEN", 
                     "AKT1", "AKT2", "AKT3", "MTOR", "TSC1", "TSC2", "RPS6KA2", "RPS6KA3", "RPS6KA6", "RHEB", 
                     "SRC", "PKT2")
p53_cellcycle_genes <- c("ATM", "CHECK2", "TP53", "MDM2", "CDKN2A", "CCND1", "E2F3", "CCNB2",
                         "MYC")
fat_cadherins_genes <- c("FAT1", "FAT2", "FAT3", "FAT4")

## genes mutated in ccRCC
genes_mutated_in_ccRCC <- c(ubiquitin_proteasome_genes,
                            swisnf_genes,
                            other_epigeneticregulator_genes,
                            pi3k_mtor_genes,
                            p53_cellcycle_genes,
                            fat_cadherins_genes)
genes_mutated_in_ccRCC <- unique(genes_mutated_in_ccRCC)


# genes for expression analysis -------------------------------------------
## RAS-pi3k-mtor
ras_pathway_genes <- c("KRAS", "HRAS", "NRAS", 
                       "BRAF", "RAF1", "MAP2K1", "MAPK1", "MAPK3")
pi3k_pathway_genes <- c("PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3", 
                        "PTEN",
                        "PDK1",
                        "AKT1", "AKT2", "AKT3",
                        "GSK3B", "BAD", "FOXO1", "FOXO3", "FOXO4", "NOS3", "CDKN1B")
mtor_pathway_genes <- c("TSC1", "TSC2", "RHEB",
                        "MTOR", "AKT1S1", "DEPTOR", "RPTOR", "MLST8", "MAPKAP1", "PRR5", "RICTOR",
                        "RB1CC1", "ULK1", "ATG13",
                        "GRB10", "TFEB", 
                        "RPS6KB1",
                        "EIF4EBP1")
## reference: https://www.wikipathways.org/index.php/Pathway:WP3844
genes_pi3k_mtor <- c(ras_pathway_genes, pi3k_pathway_genes, mtor_pathway_genes)
## RTKs
met_related_genes <- c("HGF", "MET", "AXL")
vegfr_genes <- c("FLT1", "KDR", "FLT3", "FLT4", "NRP1", "NRP2")
vegf_genes <- c("VEGFA", "VEGFB", "VEGFC", "VEGFD", "VEGFE")
## PMID: 21447729
other_cabo_related_genes <- c("KIT", "RET", "NTRK2", "TEK", "ROS1", "MERTK", "TYRO3",
                              "NTRK1", "IGF1R", "FGFR1", "SOX18", "MACC1")
## https://ascopubs.org/doi/10.1200/JCO.2006.06.3602?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub++0pubmed&
sunitinib_genes <- c("KDR", "PDGFRB", "PDGFRA", "FGFR1", "FLT3", "KIT", "CSF1R", "RET")
genes_rtk_cabo <- c(met_related_genes, vegfr_genes, vegf_genes, other_cabo_related_genes, other_cabo_related_genes, sunitinib_genes)
genes_rtk_cabo <- unique(genes_rtk_cabo)
