## set fixed parameters for Seurat
num_pcs <- 40
num_var_features <- 3000

## set prefix for feature names for human genes and mouse genes
prefix_human_feature <- "GRCh38-3.0.0.premrna-"
prefix_mouse_feature <- "mm10-premrna---------"

## SMGs
SMGs <- list()
SMGs[["CCRCC"]] <- c("VHL", "PBRM1", "SETD2", "KDM5C", "PTEN", "BAP1", "MTOR", "TP53")
