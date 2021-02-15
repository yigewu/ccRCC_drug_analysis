# Yige Wu @ WashU 2018 Jan
## shared functions for cptac2p_analysis

# set variables -----------------------------------------------------------
## base directory at Box Sync
# baseD = "~/Box/"
# setwd(baseD)
# dir2dinglab_projects <- paste0(baseD, "Ding_Lab/Projects_Current/")
# ## processed data directory
# dir_analysis_result <- paste0(dir2dinglab_projects, "RCC/ccRCC_Drug/Resources/Analysis_Results/")
dir_analysis_result <- paste0(dir_base, "/Resources/Analysis_Results/")
  
# make directory ----------------------------------------------------------
makeOutDir_katmai = function(path_script) {
  folders <- strsplit(x = path_script, split = "\\/")[[1]]
  folder_num <- which(folders == "ccRCC_drug_analysis") + 1
  dir_analysis_resultnow <- paste(strsplit(paste(folders[folder_num:length(folders)], collapse = "/"), split = "\\.")[[1]][1], sep = "/")
  dir_analysis_resultnow <- paste0(dir_analysis_result, dir_analysis_resultnow, "/")
  dir.create(dir_analysis_resultnow)
  dir_analysis_resultnow_son <- dir_analysis_resultnow
  dirs2make <- NULL
  while (!dir.exists(dir_analysis_resultnow_son)) {
    tmp <- strsplit(dir_analysis_resultnow_son, split = "\\/")[[1]]
    dir_analysis_resultnow_parent <-paste(tmp[-length(tmp)], collapse = "/")
    dir.create(dir_analysis_resultnow_parent)
    dir.create(dir_analysis_resultnow_son)
    dir.create(dir_analysis_resultnow)
    if (!dir.exists(dir_analysis_resultnow_son)) {
      dirs2make[length(dirs2make) + 1] <- dir_analysis_resultnow_son
    }
    dir_analysis_resultnow_son <- dir_analysis_resultnow_parent
  }
  
  if (length(dirs2make) > 0){
    for (i in 1:length(dirs2make)) {
      dir.create(dirs2make[i])
    }
  } 
  return(dir_analysis_resultnow)
}

makeOutDir = function() {
  folders <- strsplit(x = rstudioapi::getSourceEditorContext()$path, split = "\\/")[[1]]
  folder_num <- which(folders == "ccRCC_drug_analysis") + 1
  dir_analysis_resultnow <- paste(strsplit(paste(folders[folder_num:length(folders)], collapse = "/"), split = "\\.")[[1]][1], sep = "/")
  dir_analysis_resultnow <- paste0(dir_analysis_result, dir_analysis_resultnow, "/")
  dir.create(dir_analysis_resultnow)
  dir_analysis_resultnow_son <- dir_analysis_resultnow
  dirs2make <- NULL
  while (!dir.exists(dir_analysis_resultnow_son)) {
    tmp <- strsplit(dir_analysis_resultnow_son, split = "\\/")[[1]]
    dir_analysis_resultnow_parent <-paste(tmp[-length(tmp)], collapse = "/")
    dir.create(dir_analysis_resultnow_parent)
    dir.create(dir_analysis_resultnow_son)
    dir.create(dir_analysis_resultnow)
    if (!dir.exists(dir_analysis_resultnow_son)) {
      dirs2make[length(dirs2make) + 1] <- dir_analysis_resultnow_son
    }
    dir_analysis_resultnow_son <- dir_analysis_resultnow_parent
  }

  if (length(dirs2make) > 0){
    for (i in 1:length(dirs2make)) {
      dir.create(dirs2make[i])
    }
  }
  return(dir_analysis_resultnow)
}

# generate mutation matrix ------------------------------------------------
generate_somatic_mutation_matrix <- function(pair_tab, maf) {
  genes4mat <- unique(unlist(pair_tab))
  length(genes4mat)
  
  maf <- maf[maf$Hugo_Symbol %in% genes4mat,]
  nrow(maf)
  maf$sampID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_T", 2)[,1]
  
  mut_mat <- reshape2::dcast(data = maf, Hugo_Symbol ~ sampID, fun =  function(x) {
    variant_class <- paste0(unique(x), collapse = ",")
    return(variant_class)
  }, value.var = "Variant_Classification", drop=FALSE)
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  return(mut_mat)
}

get_somatic_mutation_vaf_matrix <- function(pair_tab, maf) {
  genes4mat <- unique(unlist(pair_tab))
  length(genes4mat)
  
  maf <- maf[maf$Hugo_Symbol %in% genes4mat,]
  nrow(maf)
  maf$sampID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
  maf$vaf <- maf$t_alt_count/(maf$t_alt_count + maf$t_ref_count)
  
  mut_mat <- reshape2::dcast(data = maf, Hugo_Symbol ~ sampID, fun =  function(x) {
    VAF <- paste0(unique(x), collapse = ",")
    return(VAF)
  }, value.var = "vaf", drop=FALSE)
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  return(mut_mat)
}

get_somatic_mutation_detailed_matrix <- function(pair_tab, maf) {
  genes4mat <- unique(unlist(pair_tab))
  length(genes4mat)
  
  maf <- maf[maf$Hugo_Symbol %in% genes4mat,]
  nrow(maf)
  maf$sampID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
  
  mut_mat <- reshape2::dcast(data = maf, Hugo_Symbol ~ sampID, fun =  function(x) {
    variant_class <- paste0(unique(x), collapse = ",")
    return(variant_class)
  }, value.var = "HGVSp_Short", drop=FALSE)
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  return(mut_mat)
}

get_somatic_mutation_aachange_vaf_matrix <- function(pair_tab, maf) {
  genes4mat <- unique(unlist(pair_tab))
  length(genes4mat)
  
  maf <- maf[maf$Hugo_Symbol %in% genes4mat & maf$Variant_Classification != "Silent",]
  nrow(maf)
  maf$sampID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
  maf$vaf <- maf$t_alt_count/(maf$t_alt_count + maf$t_ref_count)
  maf$HGVSp_sim <- gsub(x = maf$HGVSp_Short, pattern = "p\\.", replacement = "")
  
  maf$aachange_vaf <- paste0(maf$HGVSp_sim, "(", signif(x = maf$vaf, digits = 2), ")")
  mut_mat <- reshape2::dcast(data = maf, Hugo_Symbol ~ sampID, fun =  function(x) {
    aahange_vaf <- paste0(unique(x), collapse = ",")
    return(aahange_vaf)
  }, value.var = "aachange_vaf", drop=FALSE)
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  return(mut_mat)
}

