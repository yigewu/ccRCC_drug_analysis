# Yige Wu @ WashU 2020 Feb
## plot VAF changes before and after treatment

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
## input data info with treatment grouping info
## which control which be paried with the treated tumors? there are 3 controls
model_seq_status_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/write_bulk_data_status/20200522.v1/RCC_PDX_Related_Samples.Batch1_10.Data_Status.20200522.v1.tsv")
## input data info with analysis id
datainfo_batch_df <- readxl::read_xlsx(path = "./Resources/Bulk_Processed_Data/Data_Info/washU-batch9.dataInfo.202004.xlsx")
## input batch mutation calls
mut_batch_df <- fread("~/Box/Ding_Lab/Projects_Current/PDX-WashU/batch9/somaticMut/b9.somaticMut.merged.remMutNearIndel.remPDXmut.meta3.tsv", data.table = F)
mut_batch_df <- mut_batch_df %>%
  filter(grepl(pattern = "RESL", x = Tumor_Sample_Barcode))

# get sampleid.across data type to analysis id ----------------------------
## add SampleID.AcrossDataTypes to datainfo_batch_df
datainfo_batch_df <- datainfo_batch_df %>%
  filter(grepl(pattern = "RESL", x = SampleID)) %>%
  mutate(SampleID.AcrossDataType = gsub(x = SampleID, pattern = "-DNA", replacement = ""))
## add Analysis.ID to the treatment grouping info
model_seq_status_df$Analysis_ID <- mapvalues(x = model_seq_status_df$SampleID.AcrossDataType, from = datainfo_batch_df$SampleID.AcrossDataType, to = datainfo_batch_df$Analysis_ID)
## make treatment group id
datainfo_treated_group_df <- model_seq_status_df %>%
  filter(!is.na(Treatment.Days)) %>%
  filter(Analysis_ID != SampleID.AcrossDataType) %>%
  mutate(Treatment_Group = paste0(ModelID, "_", Treatment.Days, "D_Treatment"))

## set treatment group
treatment_group_tmp <- "RESL10_38D_Treatment"
for (treatment_group_tmp in unique(datainfo_treated_group_df$Treatment_Group)) {
  ## create output directory
  dir_out_sub <- paste0(dir_out, treatment_group_tmp, "/")
  dir.create(dir_out_sub)
  
  ## get analysis ids for control and treated
  analysis_id_control <- datainfo_treated_group_df$Analysis_ID[datainfo_treated_group_df$Treatment_Group == treatment_group_tmp & datainfo_treated_group_df$TumorTissue == "PDX tumor;control"]
  analysis_id_control
  analysis_id_control <- analysis_id_control[1]
  analysis_ids_treated <- datainfo_treated_group_df$Analysis_ID[datainfo_treated_group_df$Treatment_Group == treatment_group_tmp & datainfo_treated_group_df$TumorTissue == "PDX tumor;treated"]
  analysis_ids_treated
  # analysis_id_treated <- analysis_ids_treated[1]
  for (analysis_id_treated in analysis_ids_treated) {
    ## create output directory
    dir_out_sub2 <- paste0(dir_out_sub, analysis_id_treated, "/")
    dir.create(dir_out_sub2)
    
    # filter maf for the paired sample --------------------------------------------
    maf_tab <- mut_batch_df %>%
      filter(Tumor_Sample_Barcode %in% c(analysis_id_treated, analysis_id_control))
    
    # get VAF for mutations in this pair ----------------------------------------------------------
    genes_allmut <- unique(maf_tab$Hugo_Symbol)
    vaf_mat <- get_somatic_mutation_vaf_matrix(pair_tab = genes_allmut, maf = maf_tab)
    vaf_mat[vaf_mat == ""] <- 0
    mut_mat <- get_somatic_mutation_detailed_matrix(pair_tab = genes_allmut, maf = maf_tab)
    mut_mat <- mut_mat %>%
      rename(control_mutation = analysis_id_control) %>%
      rename(treated_mutation = analysis_id_treated) %>%
      mutate(mutation = ifelse(control_mutation == "", treated_mutation, control_mutation))
    vaf_mat$mutation <- mut_mat$mutation
    file2write <- paste0(dir_out_sub2, "Somatic_Mutation_VAF_Changes.Cabo_Treated_vs_Control.", run_id, ".tsv")
    write.table(x = vaf_mat, file = file2write, quote = F, sep = "\t", row.names = F)
    
    # plot scatterplot --------------------------------------------------------
    plot_data <- vaf_mat %>%
      rename(control_vaf = analysis_id_control) %>%
      rename(treated_vaf = analysis_id_treated)
    ## make sure the numbers are in numeric format
    plot_data$control_vaf <- as.numeric(as.vector(plot_data$control_vaf))
    plot_data$treated_vaf <- as.numeric(as.vector(plot_data$treated_vaf))
    ## add text to show in the color legend
    plot_data <- plot_data  %>%
      mutate(color = ifelse(treated_vaf > (control_vaf + 0.1), "treated_vaf > control_vaf + 10%", 
                            ifelse(treated_vaf < (control_vaf - 0.1), "treated_vaf < control_vaf - 10%", "|treated_vaf - control_vaf| < 10%")))
    
    ## plot mutations in SMGs
    ### filter the data frame to get the ones to highlight
    text_data <- plot_data %>%
      filter((Hugo_Symbol %in% ccRCC_SMGs))
    
    p <- ggplot()
    p <- p + geom_point(data = plot_data, mapping = aes(x = control_vaf, y = treated_vaf))
    p <- p + geom_abline(slope = 1, linetype = 2)
    p <- p + geom_text_repel(data = text_data, mapping = aes(x = control_vaf, y = treated_vaf, label = Hugo_Symbol, color = color))
    p <- p + scale_color_manual(values = c("treated_vaf > control_vaf + 10%" = "red", "treated_vaf < control_vaf - 10%" = "blue", "|treated_vaf - control_vaf| < 10%" = "black"))
    p <- p + theme(legend.position = "top")
    file2write <- paste0(dir_out_sub2, "Somatic_Mutation_VAF_Changes.", analysis_id_treated, "_vs_Control.", "SMGs", ".png")
    png(file = file2write, width = 800, height = 800, res = 150)
    print(p)
    dev.off()
    
    ## plot mutations in with VAF Increase
    ### filter the data frame to get the ones to highlight
    text_data <- plot_data %>%
      filter(color == "treated_vaf > control_vaf + 10%")
    
    p <- ggplot()
    p <- p + geom_point(data = plot_data, mapping = aes(x = control_vaf, y = treated_vaf))
    p <- p + geom_abline(slope = 1, linetype = 2)
    p <- p + geom_text_repel(data = text_data, mapping = aes(x = control_vaf, y = treated_vaf, label = Hugo_Symbol, color = color))
    p <- p + scale_color_manual(values = c("treated_vaf > control_vaf + 10%" = "red", "treated_vaf < control_vaf - 10%" = "blue", "|treated_vaf - control_vaf| < 10%" = "black"))
    p <- p + theme(legend.position = "top")
    file2write <- paste0(dir_out_sub2, "Somatic_Mutation_VAF_Changes.", analysis_id_treated, "_vs_Control.", "Increasted_VAF", ".png")
    png(file = file2write, width = 800, height = 800, res = 150)
    print(p)
    dev.off()
    
    ## plot mutations in with VAF Decrease
    ### filter the data frame to get the ones to highlight
    text_data <- plot_data %>%
      filter(color == "treated_vaf < control_vaf - 10%")
    
    p <- ggplot()
    p <- p + geom_point(data = plot_data, mapping = aes(x = control_vaf, y = treated_vaf))
    p <- p + geom_abline(slope = 1, linetype = 2)
    p <- p + geom_text_repel(data = text_data, mapping = aes(x = control_vaf, y = treated_vaf, label = Hugo_Symbol, color = color))
    p <- p + scale_color_manual(values = c("treated_vaf > control_vaf + 10%" = "red", "treated_vaf < control_vaf - 10%" = "blue", "|treated_vaf - control_vaf| < 10%" = "black"))
    p <- p + theme(legend.position = "top")
    file2write <- paste0(dir_out_sub2, "Somatic_Mutation_VAF_Changes.", analysis_id_treated, "_vs_Control.", "Decreased_VAF", ".png")
    png(file = file2write, width = 800, height = 800, res = 150)
    print(p)
    dev.off()
  }
}
