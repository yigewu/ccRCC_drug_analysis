# Yige Wu @ WashU 2019 Aug
## plot a heatmap with genomics data and proteomics data and kinase-substrate score status for given pairs


###########################################
######## Source
###########################################
# source ------------------------------------------------------------------
setwd(dir = "/Users/yigewu/Box/")
source(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_drug/pdx_r_analysis/pdx_shared.R")

# set variables -----------------------------------------------------------
## plotting paramters
cap <- 3
breaks = seq(-(cap),cap, by=0.2)
## add color palette
color.palette <- colorRampPalette(rev(brewer.pal(10, "RdBu")))(length(breaks))
color.palette.OrRd <- colorRampPalette(brewer.pal(10, "OrRd"))(length(breaks))
color.palette.Bu <- colorRampPalette(rev(brewer.pal(10, "Blues")))(length(breaks))[1:15]

# plot VHL Expressions for CCRCC-----------------------------------------------
mut_genes <- c(SMGs[["CCRCC"]])
cna_genes <- c("CDKN2A", "VHL")
rna_genes <- c("HIF1A","EPAS1", 
               "VEGFA", ## VEGFA controlled by HIFs
               "FGF2", ## mostly FGF2 is reported
               "FGFR1",
               "MET", "AXL",
               "FLT4", "FLT1", "KDR", "FLT3", 
               "PDGFA","PDGFB", ## PDGFB controlled by /HIFs
               "PDGFRA", "PDGFRB",
               "MTOR", "PTEN", "VHL", "PBRM1", "BAP1", "SETD2",
               "PDCD1", "CD274", "PDCD1LG2", "CTLA4",
               "CDK4", "CDK6", "CDK7")

pro_genes <- c("")
phog_genes <- c("")
pho_genes_uniq <- c( "")
pho_genes <- pho_genes_uniq
rsds <- c("")
row_order <- c(paste0(rna_genes, "_RNA"))
fig_width <- 10
fig_height <- 3
nonNA_cutoff <- 0
version_tmp <- 3
if_cluster_row_tmp <- F
if_cluster_col_tmp <- F
if_col_name_tmp <- T
is_cna_quantitative <- T

# bussiness ------------------------------------------------------------------
geneA <- paste(head(unique(c(mut_genes, cna_genes)), 5), collapse = "_")
geneB <- paste(head(unique(c(rna_genes, pro_genes, pho_genes)), 5), collapse = "_")
Expression <- paste0(rsds, collapse = "_")

# for (cancer in c("UCEC", "BRCA", "CCRCC", "CO", "OV")) {
for (cancer in c("CCRCC")) {
  subdir1 <- paste0(makeOutDir(), cancer, "/")
  dir.create(subdir1)
  
  fn <- paste0(subdir1, paste(geneA, geneB, Expression, sep = "_"), "_", cancer, ".pdf")
  
  if (!file.exists(fn)) {
    ann_colors <- list()
    
    # input data first because different for each cancer type --------------------------------------------------------------
    ## input mutation matrix
    maf <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_PDX/Analysis/summary_results/wxs.somaticMut/somaticMut_merged.20190714/u54.wxs_somaticMut_merged.maf.rc.caller", data.table = F)
    maf <- maf %>%
      filter(Tumor_Sample_Barcode %in% c("PDX_WUR_014_T", "PDX_WUR_016_T", "PDX_WUR_065_3548_T", "PDX_WUR_065_3549_T")) %>%
      filter(Hugo_Symbol %in% mut_genes)
    mut_mat <- generate_somatic_mutation_matrix(pair_tab = mut_genes, maf = maf)
    
    ## load RNA
    rna_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_PDX/Analysis/summary_results/rna.geneExp/u54.trans2gene/b1-b6.trans2gene/kallisto.tpm.gene_level.tsv", data.table = F)
    rna_tab <- rna_tab[rna_tab$Gene %in% rna_genes, c("Gene", "TWDE-WUR-014-TRESL5B_4495m", "TWDE-WUR-016-TRESL3C_3317m", "TWDE-WUR-065-TRESL4C_3549m", "TWDE-WUR-065-TRESL4C_3548m")]

    # load CNA ----------------------------------------------------------------
    cnv_file_names <- list.files(path = "./Ding_Lab/Projects_Current/RCC/ccRCC_PDX/Analysis/summary_results/wxs.cnv.ccRCC/orig.result.20190714/")
    cnv_file_names <- cnv_file_names[grepl(pattern = "tsv", x = cnv_file_names)]
    cnv_file_names
    cna_tab <- NULL
    for (cnv_file_name in cnv_file_names) {
      cnv_tab_tmp <- fread(input = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_PDX/Analysis/summary_results/wxs.cnv.ccRCC/orig.result.20190714/", cnv_file_name))
      cnv_tab_tmp$New_ID <- str_split_fixed(string = cnv_file_name, pattern = ".T.segment_gene.chr.tsv", n = 2)[,1]
      cna_tab <- rbind(cna_tab, cnv_tab_tmp %>%
                         dplyr::filter(gene %in% cna_genes))
    }
    cna_tab <- cna_tab %>%
      select(gene, log2, New_ID)
    
    # make the annotation columns for each sample -----------------------------
    col_anno <- data.frame(New_ID = c("PDX_WUR_014", "PDX_WUR_016", "PDX_WUR_065_3548", "PDX_WUR_065_3549"))
    
    if (nrow(mut_mat) > 0){
      mut_mat.m <- melt(mut_mat, id.vars = "Hugo_Symbol")
      mut_mat.m %>% head()
      mut_mat.m <- data.frame(mut_mat.m)
      colnames(mut_mat.m) <- c("Gene", "New_ID", "variant_class")
      
      ## distinguish by missense and truncation
      mut_mat.m$variant_class[is.na(mut_mat.m$variant_class)] <- ""
      mut_mat.m$variant_class_sim <- "other_mutation"
      mut_mat.m$variant_class_sim[mut_mat.m$variant_class == ""] <- "wild_type"
      mut_mat.m$variant_class_sim[mut_mat.m$variant_class  == "Silent"] <- "silent"
      mut_mat.m$variant_class_sim[grepl(x = mut_mat.m$variant_class, pattern = "Missense_Mutation")] <- "missense"
      mut_mat.m$variant_class_sim[grepl(x = mut_mat.m$variant_class, pattern = "Nonsense_Mutation|Frame_Shift_Ins|Frame_Shift_Del")] <- "truncation"
      mut_mat.m$variant_class_sim[sapply(X = mut_mat.m$variant_class, FUN = function(v) (grepl(pattern = "Nonsense_Mutation|Frame_Shift_Ins|Frame_Shift_Del", x = v) & grepl(pattern = "Missense_Mutation", x = v)))] <- "missense&truncation"
      
      for (gene in unique(mut_mat.m$Gene[mut_mat.m$variant_class_sim != "wild_type"])) {
        mut_mat2merge <- mut_mat.m[mut_mat.m$Gene == gene, c("New_ID", "variant_class_sim")]
        colnames(mut_mat2merge) <- c("New_ID", paste0("mutation.", gene))
        col_anno <- merge(col_anno, mut_mat2merge, by = c("New_ID"), all.x = T)
      }
      
      for (gene in mut_genes) {
        if (paste0("mutation.", gene) %in% colnames(col_anno)) {
          col_anno <- col_anno[order(col_anno[, paste0("mutation.", gene)], decreasing = T),]
          ann_colors[[paste0("mutation.", gene)]] <- c(missense = "#E41A1C", truncation = "#377EB8", wild_type = "white", "missense&truncation" = "#6A3D9A", other_mutation = "#FF7F00", silent = "#33A02C")
        }
      }
      
    } else {
      print("no mutation!")
    }
    
    ## CNA needs to show both geneA and geneB; levels: amplification, deletion, neutral
    if (nrow(cna_tab) > 0) {
      cna_tab.m <- cna_tab
      if (is_cna_quantitative == T) {
        for (gene in intersect(cna_genes, cna_tab.m$gene)) {
          cna_mat2merge <- cna_tab.m[cna_tab.m$gene == gene, c("New_ID", "log2")]
          colnames(cna_mat2merge) <- c("New_ID", paste0("log2.", gene))
          col_anno <- merge(col_anno, cna_mat2merge, by = c("New_ID"), all.x = T)
        }
        for (gene in cna_genes) {
          if (paste0("log2.", gene) %in% colnames(col_anno)) {
            col_anno <- col_anno[order(col_anno[, paste0("log2.", gene)], decreasing = T),]
            ann_colors[[paste0("log2.", gene)]] <-  color.palette.Bu
          }
        }
        
      } else {
        for (gene in intersect(cna_genes, unique(cna_tab.m$gene[cna_tab.m$CNA != "neutral"]))) {
          cna_mat2merge <- cna_tab.m[cna_tab.m$gene == gene, c("New_ID", "CNA")]
          colnames(cna_mat2merge) <- c("New_ID", paste0("CNA.", gene))
          col_anno <- merge(col_anno, cna_mat2merge, by = c("New_ID"), all.x = T)
        }
        
        for (gene in cna_genes) {
          if (paste0("CNA.", gene) %in% colnames(col_anno)) {
            col_anno <- col_anno[order(col_anno[, paste0("CNA.", gene)], decreasing = T),]
            ann_colors[[paste0("CNA.", gene)]] <-  c(amplification = "#E41A1C", deletion = "#377EB8", "neutral" = "grey")
          }
        }
      }
    } else {
      print("no CNA!")
    }
    
    
    col_anno %>% head()
    rownames(col_anno) <- col_anno$New_ID
    col_anno$New_ID <- NULL
    
    # make the matrix of values showing in heatmap ----------------------------
    sample_ID_map <- readxl::read_excel(path = "./Ding_Lab/Projects_Current/RCC/ccRCC_PDX/Analysis/summary_data_info/all_sample_info.v20190701.xlsx")
    sample_ID_map <- data.frame(sample_ID_map)
    sup_tab_can <- NULL
    
    if (nrow(rna_tab) > 0) {
      rna_tab.m <- melt(rna_tab, id.vars = "Gene")
      colnames(rna_tab.m) <- c("Gene", "variable", "TPM")
      ## shitty one letter off from the MGI_ID.v2
      rna_tab.m$MGI_ID.v2 <- paste0(str_split_fixed(string = rna_tab.m$variable, pattern = "TRESL", n = 2)[,1], "DRESL", str_split_fixed(string = rna_tab.m$variable, pattern = "TRESL", n = 2)[,2])
      rna_tab.m <- rna_tab.m %>%
        mutate(exp_value = log2(TPM+1)) %>%
        mutate(Expression = "RNA")
      rna_tab.m <- merge(rna_tab.m, sample_ID_map %>%
                           select(MGI_ID.v2, New_ID.v1), by = c("MGI_ID.v2"), all.x = T)
      
      rna_tab.m <- rna_tab.m %>%
        mutate(New_ID = New_ID.v1) %>%
        filter(New_ID.v1 %in% rownames(col_anno)) %>%
        unique()
      sup_tab_can <- rbind(sup_tab_can, rna_tab.m[,c("Gene", "Expression", "New_ID", "exp_value")])
    }
    
    sup_tab_can$id_row <- paste0(sup_tab_can$Gene, "_", sup_tab_can$Expression)
    sup_tab_can$exp_value <- as.numeric(as.vector(sup_tab_can$exp_value))
    sup_tab_can <- unique(sup_tab_can)
    
    ## make the matrix for the heatmap body
    df_value <- dcast(data = sup_tab_can, New_ID ~ id_row, value.var = "exp_value")
    
    df_value %>% head()
    mat_value <- as.matrix(df_value[,-1])
    rownames(mat_value) <- df_value$New_ID
    mat_value %>% head()
    
    ## order the matrix column
    
    ## order the matrix rows
    # if (length(row_order) > 1) {
    #   mat_value <- mat_value[intersect(row_order, rownames(mat_value)),]
    #   mat_value <- mat_value[rowSums(!is.na(mat_value)) >= nonNA_cutoff, ]
    # } else {
    #   mat_value <- matrix(data = mat_value, nrow = 1, dimnames = list(row_order, names(mat_value)))
    # }
    
    
    fn <- paste0(makeOutDir(), paste(geneB, Expression, sep = "_"), ".", format(Sys.Date(), "%Y%m%d") , ".v", version_tmp, ".pdf")
    my_heatmap <- pheatmap(mat_value, 
                           color = color.palette.OrRd,
                           annotation_row = col_anno,
                           annotation_colors = ann_colors,
                           na_col = "white",
                           show_colnames = if_col_name_tmp,
                           cluster_rows=if_cluster_row_tmp, 
                           cluster_cols=if_cluster_col_tmp)
    save_pheatmap_pdf(x = my_heatmap, 
                      filename = fn, 
                      width = fig_width, height = fig_height)
    
  }
}


