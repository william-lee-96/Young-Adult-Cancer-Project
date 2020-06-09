##### ImmuneFeat_YoungOnset_Assoc.R #####
# William Lee @ 2019
# conduct association between immune features and young onset cancer (<50 or age as continuous variable)
# load "readxl" and "tidyverse"
# set working directory for dev and debug

bdir = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/immune_features"
setwd(bdir)

# read in statistical functions for glm wrapper
source("../stat_functions.R")
system("mkdir out")

library(readxl)
library(tidyverse)

# reads excel spreadsheet into R dataframe
immune_feat_xlfile = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/immune_features/1-s2.0-S1074761318301213-mmc2.xlsx"
immune_feat_df = read_xlsx(immune_feat_xlfile, sheet = NULL, range = NULL, col_names = TRUE,
                          col_types = NULL, na = "NA", trim_ws = TRUE, skip = 0,
                          progress = readxl_progress(), .name_repair = "unique")
unique_immune_samples = unique(substr(immune_feat_df$"TCGA Participant Barcode",1,12))

colnames(immune_feat_df)[1] <- 'bcr_patient_barcode'
colnames(immune_feat_df)[2] <- 'TCGA_Study'
colnames(immune_feat_df)[5:length(colnames(immune_feat_df))] <- gsub('[[:space:]]', '_', colnames(immune_feat_df)[5:length(colnames(immune_feat_df))])
colnames(immune_feat_df)[5:length(colnames(immune_feat_df))] <- gsub('[[:punct:]]', '_', colnames(immune_feat_df)[5:length(colnames(immune_feat_df))])

immune_features = colnames(immune_feat_df)[5:length(colnames(immune_feat_df))]

##### clinical files #####
# clinical file #
clin_complete_f = "~/Box/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/GDC_Data_2018/clinical_PANCAN_patient_with_followup.tsv"
clin_complete = read.table(header=T, quote = "", sep="\t", fill =T, file = clin_complete_f, stringsAsFactors=FALSE)

# Principal Component (PC) file #
PCs_f = "~/Box/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/Germline_Huang_Cell2018/2017-04-24_GDAN_AIM_PCA_ethnicity_assigned_WashU.tsv"
PCs = read.table(header=T, quote = "", sep="\t", fill =T, file = PCs_f, stringsAsFactors=FALSE)

clin_brief = clin_complete[,c("bcr_patient_barcode","acronym","age_at_initial_pathologic_diagnosis","gender")]
colnames(PCs)[1] = "bcr_patient_barcode"
PCs = PCs[!duplicated(PCs$bcr_patient_barcode),]
clin_merge = merge(clin_brief,PCs, by= "bcr_patient_barcode",all.x=T,all.y=F)

# subtype file #
subtype_f = "~/Box/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/GDC_Data_2018/Subtype\ Assignments.xlsx"
subtype = data.frame(readxl::read_xlsx(subtype_f))
colnames(subtype)[1] = "bcr_patient_barcode"
clin_merge_subtype = merge(clin_merge,subtype, by= "bcr_patient_barcode")

# create binary age variable; typically <=50 considered as young onset cancer
clin_merge_subtype$age_binary = clin_merge_subtype$age_at_initial_pathologic_diagnosis <= 50

# only include the ones with available molecular data
clin_merge_subtype_avail = clin_merge_subtype[clin_merge_subtype$bcr_patient_barcode %in% unique_immune_samples,]
clin_merge_subtype_avail_complete = clin_merge_subtype_avail[complete.cases(clin_merge_subtype_avail[,2:4]) & clin_merge_subtype_avail$gender != "",]

# quick check of merged data
cat("TCGA cases count available for this analysis\n")
nrow(clin_merge_subtype_avail_complete)
cat("Cancer type distribution of TCGA cases with young onset cancer\n")
table(clin_merge_subtype_avail_complete$age_at_initial_pathologic_diagnosis<=50,clin_merge_subtype_avail_complete$acronym)

##### statistical testing: immune feature analysis within each cancer type #####
# immune features ~ onset age + covariates ("gender","Subtype","PC1","PC2")

clin_merge_subtype_avail_complete$cancer = NULL
clin_merge_subtype_avail_complete$ethnicity = NULL
clin_merge_subtype_avail_complete$washu_assigned_ethnicity = NULL
clin_merge_subtype_avail_complete$Sample = NULL
clin_merge_subtype_avail_complete[ clin_merge_subtype_avail_complete == "NA" ] <- NA
clin_merge_subtype_avail_complete[ clin_merge_subtype_avail_complete == "[Not Available]"] <- NA

# initiate empty list and index = 1 to store results
immune_results_list = as.list(NULL)
immune_results_list_binary = as.list(NULL)
i=1

# conduct tests by each cancer (denoted as acronym)
for (cancer in unique(clin_merge_subtype_avail_complete$acronym)) {
  
  if (cancer == "BRCA" | cancer == "CESC" | cancer == "COAD" | cancer == "HNSC" |
      cancer == "KIRC" | cancer == "KIRP" | cancer == "LGG"  | cancer == "LIHC" |
      cancer == "OV"   | cancer == "PCPG" | cancer == "SARC" | cancer == "SKCM" | 
      cancer == "THCA" | cancer == "UCEC") {
  
  clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete[clin_merge_subtype_avail_complete$acronym==cancer,]
  clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete_c[complete.cases(clin_merge_subtype_avail_complete_c),]
  
  # conduct the test by immune feature
  for (immune_feat in immune_features) {
  
   # if (immune_feat != "OS" & immune_feat != "PFI" & immune_feat != "BCR_Evenness" & immune_feat != "BCR_Shannon" & immune_feat != "BCR_Richness" &
   #     immune_feat != "TCR_Evenness" & immune_feat != "TCR_Shannon" & immune_feat != "TCR_Richness") {
    
   # if (immune_feat == "Th1_Cells" | immune_feat == "Th2_Cells" | immune_feat == "Th17_Cells") {
     
   # if (immune_feat == "Lymphocyte_Infiltration_Signature_Score" | immune_feat == "Wound_Healing" | immune_feat == "Proliferation" |
   #     immune_feat == "Macrophage_Regulation" | immune_feat == "IFN_gamma_Response" | immune_feat == "TGF_beta_Response") {
    
     if (immune_feat == "Macrophages" | immune_feat == "Macrophages_M1" | immune_feat == "Macrophages_M2" |
         immune_feat == "Dendritic_Cells" | immune_feat == "B_Cells_Naive" | immune_feat == "T_Cells_CD4_Naive" |
         immune_feat == "T_Cells_CD8" | immune_feat == "NK_Cells_Activated" | immune_feat == "NK_Cells_Resting" |
         immune_feat == "T_Cells_Regulatory_Tregs") {
    
    # if (immune_feat == "Indel_Neoantigens" | immune_feat == "SNV_Neoantigens") {
    #   if (cancer == "SKCM" & immune_feat == "SNV_Neoantigens") {next}
    #   if (cancer == "SKCM" & immune_feat == "Indel_Neoantigens") {next}
    #   if (cancer == "OV" & immune_feat == "Indel_Neoantigens") {next}
    
    # if (immune_feat == "Silent_Mutation_Rate" | immune_feat == "Nonsilent_Mutation_Rate") {
       
      immune_feat_df_current = select(immune_feat_df, bcr_patient_barcode, immune_feat)
      immune_feat_df_current = na.omit(immune_feat_df_current)

      clin_merge_subtype_avail_complete_c_merge_if = merge(clin_merge_subtype_avail_complete_c,immune_feat_df_current, by= "bcr_patient_barcode", all.x = T, all.y = F)
      clin_merge_subtype_avail_complete_c_merge_if = clin_merge_subtype_avail_complete_c_merge_if[complete.cases(clin_merge_subtype_avail_complete_c_merge_if[29]),]

      if (dim(clin_merge_subtype_avail_complete_c_merge_if)[1] == 0) {next}

      if (all(clin_merge_subtype_avail_complete_c_merge_if$age_binary) | 
          any(clin_merge_subtype_avail_complete_c_merge_if$age_binary) == FALSE) {next}
      
      else {
      ## model onset age as a linear variable: TODO: add in gender as a covariate 
      model_results = run_glm(data = clin_merge_subtype_avail_complete_c_merge_if, yi = immune_feat, xi = "age_at_initial_pathologic_diagnosis", ytype = "Continuous", covi = c("SUBTYPE","gender","PC1","PC2"))
      cancer_stat = data.frame(cbind(cancer, immune_feat, model_results))
      immune_results_list[[i]] = cancer_stat
      # 
      # model onset age as a binary variable
      model_results = run_glm(data = clin_merge_subtype_avail_complete_c_merge_if, yi = immune_feat, xi = "age_binary", ytype = "Continuous", covi = c("SUBTYPE","gender","PC1","PC2"))
      cancer_stat = data.frame(cbind(cancer, immune_feat, model_results))
      immune_results_list_binary[[i]] = cancer_stat
      }
    }
    
    # increment index
    i = i + 1
  }
 }
}


##### Compile and store results #####
# compile linear result
# tt = do.call(rbind,immune_results_list)
# colnames(tt) = c("cancer","immune_feature","y","y_type","x","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
#                  "F", "Pr_F", "coefficient","covariates");
# # multiple-testing correction and sort by significance
# tt$FDR = p.adjust(tt[,"Pr_F"], method="fdr") 
# tt=tt[order(tt$Pr_F, decreasing=FALSE),]
# tn = "out/immune_features_onsetAge_assoc.txt"
# write.table(tt, quote=F, sep="\t", file = tn, row.names = F)

# compile binary result
tt = do.call(rbind,immune_results_list_binary)
colnames(tt) = c("cancer","immune_feature","y","y_type","x","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
                 "F", "Pr_F","coefficient","covariates");
# multiple-testing correction and sort by significance
tt$FDR = p.adjust(tt[,"Pr_F"], method="fdr") 
tt=tt[order(tt$Pr_F, decreasing=FALSE),]
#tn = "out/immune_features_onsetAgeBinary_assoc.txt"
#tn = "out/Th_Cells_onsetAgeBinary_assoc.txt"
#tn = "out/IGS_onsetAgeBinary_assoc.txt"
tn = "out/IIC_onsetAgeBinary_assoc.txt"
#tn = "out/Neoantigens_onsetAgeBinary_assoc.txt"
#tn = "out/Mut_Rate_onsetAgeBinary_assoc.txt"
write.table(tt, quote=F, sep="\t", file = tn, row.names = F)

#
#
#
#
#
#
#
#
#
#

#-----------Code that might be useful at some point-----------#

# if (cancer == "SKCM" & immune_feat == "CTA_Score") {next}
# if (cancer == "CHOL" & immune_feat == "BCR_Evenness") {next}
# if (cancer == "CHOL" & immune_feat == "BCR_Shannon") {next}
# if (cancer == "CHOL" & immune_feat == "BCR_Richness") {next}
# if (cancer == "LGG" & immune_feat == "BCR_Evenness") {next}
# if (cancer == "KICH" & immune_feat == "BCR_Evenness") {next}
# if (cancer == "KICH" & immune_feat == "BCR_Shannon") {next}
# if (cancer == "KICH" & immune_feat == "BCR_Richness") {next}
# if (cancer == "ACC" & immune_feat == "BCR_Evenness") {next}
# if (cancer == "ACC" & immune_feat == "BCR_Shannon") {next}
# if (cancer == "ACC" & immune_feat == "BCR_Richness") {next}
# if (cancer == "ACC" & immune_feat == "TCR_Evenness") {next}
# if (cancer == "LGG" & immune_feat == "BCR_Evenness") {next}
# if (cancer == "PCPG" & immune_feat == "BCR_Evenness") {next}
# if (cancer == "UVM" & immune_feat == "BCR_Evenness") {next}
# if (cancer == "UVM" & immune_feat == "BCR_Shannon") {next}
# if (cancer == "UVM" & immune_feat == "BCR_Richness") {next}

#violin_brca_Th17_data = select(clin_merge_subtype_avail_complete_c_merge_if, acronym, age_binary, Th17_Cells)
#violin_sarc_Th17_data = select(clin_merge_subtype_avail_complete_c_merge_if, acronym, age_binary, Th17_Cells)
#violin_ucec_Th17_data = select(clin_merge_subtype_avail_complete_c_merge_if, acronym, age_binary, Th17_Cells)

#violin_brca_Th1_data = select(clin_merge_subtype_avail_complete_c_merge_if, acronym, age_binary, Th1_Cells)
#violin_sarc_Th1_data = select(clin_merge_subtype_avail_complete_c_merge_if, acronym, age_binary, Th1_Cells)
#violin_ucec_Th1_data = select(clin_merge_subtype_avail_complete_c_merge_if, acronym, age_binary, Th1_Cells)

#violin_brca_snv_data = select(clin_merge_subtype_avail_complete_c_merge_if, acronym, age_binary, SNV_Neoantigens)
#violin_sarc_snv_data = select(clin_merge_subtype_avail_complete_c_merge_if, acronym, age_binary, SNV_Neoantigens)
#violin_ucec_snv_data = select(clin_merge_subtype_avail_complete_c_merge_if, acronym, age_binary, SNV_Neoantigens)

# separate out male/female, young vs. not.
# focus on the plots that look good/significant (do them all if you get the chance?)
# color ucec by msi status (see subtype column)
# check 15 cancers to make sure they have subtype information (NAs are automatically dropped by glm analysis)

#   ethnicity = clin_merge_subtype_avail_complete_c$washu_assigned_ethnicity
#   age_binary = clin_merge_subtype_avail_complete_c$age_binary
#   print(cancer)
#   print(table(ethnicity, age_binary))
#   #for (ethnicity in washu_ethnicity) {
#   }
# }

###########################################################################above: demographic table code

#clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete[clin_merge_subtype_avail_complete$acronym==cancer,]
#clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete_c[!is.na(clin_merge_subtype_avail_complete_c$age_binary),]
#clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete_c[!is.na(clin_merge_subtype_avail_complete_c$washu_assigned_ethnicity),]

#print(cancer)
#print(table(clin_merge_subtype_avail_complete_c$age_binary))

# it worked!