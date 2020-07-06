##### YoungAdult_Amplification_Clin_Assoc.R #####
# William Lee @ January 2020
# Updated June 2020

bdir = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/methylation_fusion"
setwd(bdir)

library(readxl)
library(tidyverse)

onc_sig_file = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/methylation_fusion/oncosigfilefixed.txt"
oncosig_df = read.table(sep="\t",header=T, file=onc_sig_file , stringsAsFactors=FALSE)

oncosig_df$SAMPLE_BARCODE <- substr(oncosig_df$SAMPLE_BARCODE, 1,12) # make sample barcode into patient barcode

names(oncosig_df)[1] <- "bcr_patient_barcode"
oncosig_df$bcr_patient_barcode= as.character(oncosig_df$bcr_patient_barcode)
unique_oncosig_samples <- unique(oncosig_df$bcr_patient_barcode) # 9125

amp_df <- cbind(oncosig_df[,1], select(oncosig_df, starts_with("AMP"))) 
names(amp_df)[1]='bcr_patient_barcode'

## clinical files
clin_complete_f = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/methylation_fusion/clinical_PANCAN_patient_with_followup.tsv" #for analysis locally
clin_complete = read.table(header=T, quote = "", sep="\t", fill =T, file = clin_complete_f, stringsAsFactors=FALSE)

## PCA files
PCs_f = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/methylation_fusion/2017-04-24_GDAN_AIM_PCA_ethnicity_assigned_WashU.tsv" #for analysis locally
PCs = read.table(header=T, quote = "", sep="\t", fill =T, file = PCs_f, stringsAsFactors=FALSE)

clin_brief = clin_complete[,c("bcr_patient_barcode","acronym","age_at_initial_pathologic_diagnosis","gender")]
colnames(PCs)[1] = "bcr_patient_barcode"
PCs = PCs[!duplicated(PCs$bcr_patient_barcode),]
clin_merge = merge(clin_brief,PCs, by= "bcr_patient_barcode",all.x=T,all.y=F) # 10956 x 28

subtype_f = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/methylation_fusion/Subtype\ Assignments.xlsx" #for analysis locally
subtype = data.frame(readxl::read_xlsx(subtype_f))
colnames(subtype)[1] = "bcr_patient_barcode"
clin_merge_subtype = merge(clin_merge,subtype, by= "bcr_patient_barcode") # 8943 x 31

clin_merge_subtype$age_binary = clin_merge_subtype$age_at_initial_pathologic_diagnosis <= 50

# only include the ones with available data in oncosig file
clin_merge_subtype_avail = clin_merge_subtype[clin_merge_subtype$bcr_patient_barcode %in% unique_oncosig_samples,] # 8943 x 32 
clin_merge_subtype_avail_complete = clin_merge_subtype_avail[complete.cases(clin_merge_subtype_avail[,2:4]) & clin_merge_subtype_avail$gender != "",] # 8943 x 32

amp_events = names(amp_df)[2:length(names(amp_df))]

clin_merge_subtype_avail_complete$cancer = NULL; clin_merge_subtype_avail_complete$ethnicity = NULL
clin_merge_subtype_avail_complete$Sample = NULL; clin_merge_subtype_avail_complete$washu_assigned_ethnicity = NULL

clin_merge_subtype_avail_complete[ clin_merge_subtype_avail_complete == "NA" ] <- NA
clin_merge_subtype_avail_complete[ clin_merge_subtype_avail_complete == "[Not Available]"] <- NA

clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="BRCA"] = "1612" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="CESC"] = "4362" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="COAD"] = "9256" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="HNSC"] = "5520" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="KIRC"] = "4467" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="KIRP"] = "4465" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="LGG"]  = "60108" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="LIHC"] = "684" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="OV"]   = "2394" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="PCPG"] = "50773" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="SARC"] = "1115" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="SKCM"] = "8923" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="THCA"] = "1781" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="UCEC"] = "1380"###

### analytic pipeline ###

# set working directory for dev and debug
bdir = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/clinical_associations" ###
setwd(bdir)

# civic_raw_data updated 06/04/2020
civic_raw_data = "01-May-2020-ClinicalEvidenceSummaries.tsv" # 3431 obs. of 41 variables
# civic_raw_data = "01-Apr-2020-ClinicalEvidenceSummaries.tsv" # 3393 obs. of 41 variables
civic_df = read.table(header=T, quote = "", sep="\t", fill=T, file = civic_raw_data, stringsAsFactors=FALSE)

### civic data preparation 
civic_df = civic_df[str_detect(civic_df$variant, "AMPLIFICATION"),]

civic_df$ampEvnt = paste(civic_df$gene, "AMP" ,sep=":")
civic_df$doid[is.na(civic_df$doid)] = "XXXXX"
### end of civic data preparation 

# cgi_raw_data updated 06/04/2020
cgi_raw_data = "cgi_biomarkers_per_variant.tsv"
cgi_df = read.table(header=T, quote = "", sep="\t", fill=T, file = cgi_raw_data, stringsAsFactors=FALSE)

### cgi data preparation
cgi_df = cgi_df[str_detect(cgi_df$Biomarker, "amplification"),]

cgi_df$Evidence.level[cgi_df$Evidence.level=="NCCN guidelines"] = "A"
cgi_df$Evidence.level[cgi_df$Evidence.level=="FDA guidelines"] = "A" 
cgi_df$Evidence.level[cgi_df$Evidence.level=="European LeukemiaNet guidelines"] = "A"
cgi_df$Evidence.level[cgi_df$Evidence.level=="NCCN/CAP guidelines"] = "A"
cgi_df$Evidence.level[cgi_df$Evidence.level=="CPIC guidelines"] = "A"

cgi_df$Evidence.level[cgi_df$Evidence.level=="Clinical trials"] = "B"
cgi_df$Evidence.level[cgi_df$Evidence.level=="Late trials"] = "B"
cgi_df$Evidence.level[cgi_df$Evidence.level=="Early trials"] = "B"
cgi_df$Evidence.level[cgi_df$Evidence.level=="Early Trials,Case Report"] = "B"

cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Glioma")] = "60108"
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "glioma")] = "60108" 
#
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Sarcoma")] = "1115"
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "sarcoma")] = "1115"
#
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Breast adenocarcinoma")] = "1612"
#
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Cervix")] = "4362"
#
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Colorectal adenocarcinoma")] = "9256"
#
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Head an neck squamous")] = "5520"
#
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Hepatic carcinoma")] = "684"
#
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Ovary")] = "2394"
#
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Cutaneous melanoma")] = "8923" 
#
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Thyroid")] = "1781" 
#
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Endometrium")] = "1380"

cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, ";")] = "XXXXX"
cgi_df$doid[is.na(cgi_df$doid)] = "XXXXX"

cgi_df$Drug = str_replace_all(cgi_df$Drug, fixed("["), "")
cgi_df$Drug = str_replace_all(cgi_df$Drug, fixed("]"), "")

cgi_df$ampEvnt = paste(cgi_df$Gene, "AMP", sep=":")
### end of cgi data preparation 

# onco_raw_data updated 06/04/2020
onco_raw_data = "oncokb_biomarker_drug_associations.tsv"
onco_df = read.table(header=T, quote = "", sep="\t", fill=T, file = onco_raw_data, stringsAsFactors=FALSE)

### onco kb data preparation
onco_df = onco_df[str_detect(onco_df$Alterations, "Amplification"),]

onco_df$Evidence.Level[onco_df$Level=="1"] = "A"
onco_df$Evidence.Level[onco_df$Level=="2"] = "A" 
onco_df$Evidence.Level[onco_df$Level=="3"] = "B"

onco_df$Evidence.Level[onco_df$Level=="R1"] = "A" 
onco_df$Evidence.Level[onco_df$Level=="R2"] = "B"

onco_df$doid[str_detect(onco_df$Tumor.Type, "Glioma")] = "60108" ###
#
onco_df$doid[str_detect(onco_df$Tumor.Type, "Breast Cancer")] = "1612" ###
#
onco_df$doid[str_detect(onco_df$Tumor.Type, "Colorectal Cancer")] = "9256" ###
#
onco_df$doid[str_detect(onco_df$Tumor.Type, "Head and Neck Squamous Cell Carcinoma")] = "5520" ###
#
onco_df$doid[str_detect(onco_df$Tumor.Type, "Melanoma")] = "8923" ###
#
onco_df$doid[str_detect(onco_df$Tumor.Type, "Liposarcoma")] = "1115" ###
onco_df$doid[str_detect(onco_df$Tumor.Type, "Ewing Sarcoma")] = "1115" ###
onco_df$doid[str_detect(onco_df$Tumor.Type, "Ewing Sarcoma of Soft Tissue")] = "1115" ###
onco_df$doid[str_detect(onco_df$Tumor.Type, "Dedifferentiated Liposarcoma")] = "1115" ###
onco_df$doid[str_detect(onco_df$Tumor.Type, "Well-Differentiated Liposarcoma")] = "1115" ###
onco_df$doid[str_detect(onco_df$Tumor.Type, "Dermatofibrosarcoma Protuberans")] = "1115" ###
#
onco_df$doid[str_detect(onco_df$Tumor.Type, "Low-Grade Serous Ovarian Cancer")] = "2394" ###
onco_df$doid[str_detect(onco_df$Tumor.Type, "Ovarian Cancer")] = "2394" ###
#
onco_df$doid[str_detect(onco_df$Tumor.Type, "Anaplastic Thyroid Cancer")] = "1781" ###
onco_df$doid[str_detect(onco_df$Tumor.Type, "Thyroid Cancer")] = "1781" ###
onco_df$doid[str_detect(onco_df$Tumor.Type, "Medullary Thyroid Cancer")] = "1781" ###
#
onco_df$doid[str_detect(onco_df$Tumor.Type, "Uterine Serous Carcinoma/Uterine Papillary Serous Carcinoma")] = "1380" ###
onco_df$doid[str_detect(onco_df$Tumor.Type, "Endometrial Cancer")] = "1380" ###

onco_df$doid[is.na(onco_df$doid)] = "XXXXX"

onco_df$ampEvnt = paste(onco_df$Gene, "AMP", sep=":")

onco_df$clin.sig[onco_df$Level=="1" | onco_df$Level=="2" | onco_df$Level=="3"] = "Responsive"
onco_df$clin.sig[onco_df$Level=="R1" | onco_df$Level=="R2"] = "Resistant"
### end of onco kb data preparation

##### the following code generates combination df fed into analytic pipeline ##### 
doid1 = civic_df$doid
doid2 = cgi_df$doid
doid3 = onco_df$doid
doid = c(doid1, doid2, doid3)

drug1 = civic_df$drugs
drug2 = cgi_df$Drug
drug3 = onco_df$Drugs
drugs = c(drug1, drug2, drug3)

ampEvnt1 = civic_df$ampEvnt
ampEvnt2 = cgi_df$ampEvnt
ampEvnt3 = onco_df$ampEvnt
ampEvnt = c(ampEvnt1, ampEvnt2, ampEvnt3)

elvl1 = civic_df$evidence_level
elvl2 = cgi_df$Evidence.level
elvl3 = onco_df$Evidence.Level
evidence_level = c(elvl1, elvl2, elvl3)

clin1 = civic_df$clinical_significance
clin2 = cgi_df$Association
clin3 = onco_df$clin.sig
clinical_significance = c(clin1, clin2, clin3)

poten_clin_act_3db = as.data.frame(cbind(doid, drugs, ampEvnt, evidence_level, clinical_significance))

######################################################################################################

unique_ampevnt_3db = unique(poten_clin_act_3db$ampEvnt)

clin_amplification_therapeutics <- data.frame(matrix(ncol = 32, nrow = 0))
colnames(clin_amplification_therapeutics) <- c("bcr_patient_barcode","acronym","age_at_initial_pathologic_diagnosis","gender",
                                        "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","SAMPLE_BARCODE",
                                        "DISEASE","SUBTYPE","age_binary","doid","amp_binary","ampEvnt","ab_drug")

# for each cancer type
for (cancer in unique(clin_merge_subtype_avail_complete$acronym)){
  
  if (cancer == "BRCA" | cancer == "CESC" | cancer == "COAD" | cancer == "HNSC" |
      cancer == "KIRC" | cancer == "KIRP" | cancer == "LGG"  | cancer == "LIHC" |
      cancer == "OV"   | cancer == "PCPG" | cancer == "SARC" | cancer == "SKCM" | 
      cancer == "THCA" | cancer == "UCEC") {
    
    clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete[clin_merge_subtype_avail_complete$acronym==cancer,]
    clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete_c[complete.cases(clin_merge_subtype_avail_complete_c),]
    
    # conduct the test by each alteration event
    for (amp in amp_events){
      
      amp_df_current = select(amp_df, bcr_patient_barcode, amp) ###
      amp_df_current = amp_df_current[amp_df_current[,amp] != 'NA',] ###
      
      clin_merge_subtype_avail_complete_c_merge_amp = merge(clin_merge_subtype_avail_complete_c, amp_df_current, by= "bcr_patient_barcode", all.x=T, all.y=F)
      clin_merge_subtype_avail_complete_c_merge_amp = clin_merge_subtype_avail_complete_c_merge_amp[clin_merge_subtype_avail_complete_c_merge_amp[,amp] != 'NA',]
      colnames(clin_merge_subtype_avail_complete_c_merge_amp)[30] <- "amp_binary"
      
      clin_merge_subtype_avail_complete_c_merge_amp$ampEvnt[clin_merge_subtype_avail_complete_c_merge_amp$amp_binary==TRUE] = paste((unlist(strsplit(amp, "[.]")))[2], "AMP", sep=":")
      clin_merge_subtype_avail_complete_c_merge_amp = clin_merge_subtype_avail_complete_c_merge_amp[complete.cases(clin_merge_subtype_avail_complete_c_merge_amp),]
      
      # skips to next loop iteration if clin_merge_subtype_avail_complete_c_merge_amp has no obs.
      if (dim(clin_merge_subtype_avail_complete_c_merge_amp)[1] == 0) {next}
      
      current_ampevnt = unique(clin_merge_subtype_avail_complete_c_merge_amp$ampEvnt)
      
      if (current_ampevnt %in% unique_ampevnt_3db) {
          
        # poten_clin_act_evnt contains only rows that have amplification event of interest 
        poten_clin_act_evnt = poten_clin_act_3db[poten_clin_act_3db$ampEvnt==current_ampevnt,]
        
        # replaces blank cells in poten_clin_act_evnt with NA (so complete.cases() will work)
        poten_clin_act_evnt[ poten_clin_act_evnt == "" ] <- NA
        
        # confirm that all rows in poten_clin_act_evnt have evidence levels and drugs
        poten_clin_act_evnt = poten_clin_act_evnt[complete.cases(poten_clin_act_evnt$evidence_level),]
        poten_clin_act_evnt = poten_clin_act_evnt[complete.cases(poten_clin_act_evnt$drugs),]
        
        # skips if poten_clin_act_evnt has no observations
        if (dim(poten_clin_act_evnt)[1] == 0) {
          clin_merge_subtype_avail_complete_c_merge_amp$ab_drug = ""
          clin_amplification_therapeutics <- rbind(clin_amplification_therapeutics, clin_merge_subtype_avail_complete_c_merge_amp)
          next
        }
        
        # instantiates empty list (elements to be appended)
        current_list = list()
          
        # initiates counter at 1
        i = 1
        
        for (elvl in poten_clin_act_evnt$evidence_level) {
          
          if (elvl == "A" & (poten_clin_act_evnt[i, "clinical_significance"] == "Sensitivity/Response" |
                             poten_clin_act_evnt[i, "clinical_significance"] == "Responsive")) {
          
          # if (elvl == "A" & (poten_clin_act_evnt[i, "clinical_significance"] == "Resistance" |
          #                    poten_clin_act_evnt[i, "clinical_significance"] == "Resistant")) {
            
            if (poten_clin_act_evnt[i, "doid"] == unique(clin_merge_subtype_avail_complete_c_merge_amp$doid)) {
              current_drug = paste(poten_clin_act_evnt[i, "drugs"], "A on-label", sep=": ")
              current_list = append(current_list, current_drug)
            }
            
            else {
              current_drug = paste(poten_clin_act_evnt[i, "drugs"], "A off-label", sep=": ")
              current_list = append(current_list, current_drug)
            }
            
            i = i + 1
            
          }
          
          else if (elvl == "B" & (poten_clin_act_evnt[i, "clinical_significance"] == "Sensitivity/Response" |
                                  poten_clin_act_evnt[i, "clinical_significance"] == "Responsive")) {
          
          # else if (elvl == "B" & (poten_clin_act_evnt[i, "clinical_significance"] == "Resistance" |
          #                         poten_clin_act_evnt[i, "clinical_significance"] == "Resistant")) {
            
            if (poten_clin_act_evnt[i, "doid"] == unique(clin_merge_subtype_avail_complete_c_merge_amp$doid)) {
              current_drug = paste(poten_clin_act_evnt[i, "drugs"], "B on-label", sep=": ")
              current_list = append(current_list, current_drug)
            }
            
            else {
              current_drug = paste(poten_clin_act_evnt[i, "drugs"], "B off-label", sep=": ")
              current_list = append(current_list, current_drug)
            }
            
            i = i + 1
            
          }
          
          else {i = i + 1}
          
        }
        
        current_list = unique(current_list)
        therapeutic_drugs = paste(current_list, collapse = " & ")
        clin_merge_subtype_avail_complete_c_merge_amp$ab_drug = therapeutic_drugs
        
        clin_amplification_therapeutics <- rbind(clin_amplification_therapeutics, clin_merge_subtype_avail_complete_c_merge_amp)
        
      }
        
      else {
        
        clin_merge_subtype_avail_complete_c_merge_amp$ab_drug = ""
        clin_amplification_therapeutics <- rbind(clin_amplification_therapeutics, clin_merge_subtype_avail_complete_c_merge_amp)
        
      }
    }
  }
}

# get dfs for actionable ampEvnt analysis later (need duplicates) #
clin_amplification_therapeutics_all = clin_amplification_therapeutics
clin_amplification_therapeutics_all = clin_amplification_therapeutics_all[
clin_amplification_therapeutics_all$acronym=="BRCA" | clin_amplification_therapeutics_all$acronym=="CESC" |
clin_amplification_therapeutics_all$acronym=="HNSC" | clin_amplification_therapeutics_all$acronym=="LGG" |
clin_amplification_therapeutics_all$acronym=="OV" | clin_amplification_therapeutics_all$acronym=="SKCM",]

clin_amplification_actionable = clin_amplification_therapeutics
clin_amplification_actionable[clin_amplification_actionable == ""] <- NA
clin_amplification_actionable = clin_amplification_actionable[!is.na(clin_amplification_actionable$ab_drug),]
clin_amplification_actionable = clin_amplification_actionable[
clin_amplification_actionable$acronym=="BRCA" | clin_amplification_actionable$acronym=="CESC" |
clin_amplification_actionable$acronym=="HNSC" | clin_amplification_actionable$acronym=="LGG" |
clin_amplification_actionable$acronym=="OV" | clin_amplification_actionable$acronym=="SKCM",]
          
# RESPONSIVE VARIANTS #
# clin_amplification_therapeutics: 3930 obs. before removing duplicates

# RESISTANT VARIANTS #
# clin_amplification_therapeutics: 3930 obs. before removing duplicates

duplicated_patient_samples = clin_amplification_therapeutics[duplicated(clin_amplification_therapeutics$bcr_patient_barcode),1]

for (duplicate in duplicated_patient_samples) {
  
  current_duplicates_df = clin_amplification_therapeutics[clin_amplification_therapeutics$bcr_patient_barcode==duplicate,]
  clin_amplification_therapeutics = clin_amplification_therapeutics[clin_amplification_therapeutics$bcr_patient_barcode!=duplicate,]
  
  drugs_as_string = paste(current_duplicates_df$ab_drug, collapse = " ")
  
  if (str_detect(drugs_as_string, "A on-label")) {
    current_duplicates_df = current_duplicates_df[str_detect(current_duplicates_df$ab_drug, "A on-label"),]
    clin_amplification_therapeutics = rbind(clin_amplification_therapeutics, current_duplicates_df[1,])
  }
  
  else if (str_detect(drugs_as_string, "A off-label")) {
    current_duplicates_df = current_duplicates_df[str_detect(current_duplicates_df$ab_drug, "A off-label"),]
    clin_amplification_therapeutics = rbind(clin_amplification_therapeutics, current_duplicates_df[1,])
  }
  
  else if (str_detect(drugs_as_string, "B on-label")) {
    current_duplicates_df = current_duplicates_df[str_detect(current_duplicates_df$ab_drug, "B on-label"),]
    clin_amplification_therapeutics = rbind(clin_amplification_therapeutics, current_duplicates_df[1,])
  }
  
  else if (str_detect(drugs_as_string, "B off-label")) {
    current_duplicates_df = current_duplicates_df[str_detect(current_duplicates_df$ab_drug, "B off-label"),]
    clin_amplification_therapeutics = rbind(clin_amplification_therapeutics, current_duplicates_df[1,])
  }
  
  else {
    clin_amplification_therapeutics = rbind(clin_amplification_therapeutics, current_duplicates_df[1,])
  }
  
}

# RESPONSIVE VARIANTS #
# clin_amplification_therapeutics: 1954 obs. after removing duplicates

# RESISTANT VARIANTS #
# clin_amplification_therapeutics: 1954 obs. after removing duplicates

### the following code produces the ggplot-ready dataframe ###

# replaces blank cells in clin_amplification_therapeutics with NA   
clin_amplification_therapeutics[clin_amplification_therapeutics == ""] <- NA

clin_amplification_therapeutics_ggplot <- data.frame(matrix(ncol = 32, nrow = 0))
colnames(clin_amplification_therapeutics_ggplot) <- colnames(clin_amplification_therapeutics)

for (cancer in unique(clin_amplification_therapeutics$acronym)){
  
  i = 1
  
  clin_amplification_therapeutics_c = clin_amplification_therapeutics[clin_amplification_therapeutics$acronym==cancer,]
  
  total_reg_cases = sum(clin_amplification_therapeutics_c$age_binary==FALSE & !is.na(clin_amplification_therapeutics_c$ab_drug))
  total_yng_cases = sum(clin_amplification_therapeutics_c$age_binary==TRUE & !is.na(clin_amplification_therapeutics_c$ab_drug))
  
  if (total_reg_cases >= 10 & total_yng_cases >= 10) {
    
    for (element in clin_amplification_therapeutics_c$ab_drug) {
      
      if (is.na(element)) {i = i + 1}
      
      else if (str_detect(element, "A on-label")) {clin_amplification_therapeutics_c[i, "ab_drug"] =  "A on-label"; i = i + 1}
      
      else if (str_detect(element, "A off-label")) {clin_amplification_therapeutics_c[i, "ab_drug"] =  "A off-label"; i = i + 1}
      
      else if (str_detect(element, "B on-label")) {clin_amplification_therapeutics_c[i, "ab_drug"] =  "B on-label"; i = i + 1}
      
      else if (str_detect(element, "B off-label")) {clin_amplification_therapeutics_c[i, "ab_drug"] =  "B off-label"; i = i + 1}
    }
    
    clin_amplification_therapeutics_ggplot <- rbind(clin_amplification_therapeutics_ggplot, clin_amplification_therapeutics_c)
    
  }
}

# RESPONSIVE VARIANTS #
# clin_amplification_therapeutics_ggplot: 1374 obs. 

# RESISTANT VARIANTS #
# clin_amplification_therapeutics_ggplot: 963 obs. 

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

clin_amplification_therapeutics_ggplot$ab_drug[is.na(clin_amplification_therapeutics_ggplot$ab_drug)] = "None"

clin_amplification_therapeutics_ggplot_rmNone = clin_amplification_therapeutics_ggplot[clin_amplification_therapeutics_ggplot$ab_drug!="None",]
clin_amplification_therapeutics_ggplot_None = clin_amplification_therapeutics_ggplot[clin_amplification_therapeutics_ggplot$ab_drug=="None",]

clin_amplification_therapeutics_ggplot$ab_drug = factor(clin_amplification_therapeutics_ggplot$ab_drug, levels=c("None", "B off-label", "B on-label", "A off-label", "A on-label"))

clin_amplification_therapeutics_ggplot$plot_age[clin_amplification_therapeutics_ggplot$age_binary==TRUE] = "<= 50"
clin_amplification_therapeutics_ggplot$plot_age[clin_amplification_therapeutics_ggplot$age_binary==FALSE] = "> 50"

clin_amplification_therapeutics_ggplot_rmNone$plot_age[clin_amplification_therapeutics_ggplot_rmNone$age_binary==TRUE] = "<= 50"
clin_amplification_therapeutics_ggplot_rmNone$plot_age[clin_amplification_therapeutics_ggplot_rmNone$age_binary==FALSE] = "> 50"

clin_amplification_therapeutics_ggplot_None$plot_age[clin_amplification_therapeutics_ggplot_None$age_binary==TRUE] = "<= 50"
clin_amplification_therapeutics_ggplot_None$plot_age[clin_amplification_therapeutics_ggplot_None$age_binary==FALSE] = "> 50"

### Selects for top amplification events ###

ampEvnt_list = list(); aeCount_list = list()

for (ae in unique(clin_amplification_therapeutics_ggplot_rmNone$ampEvnt)) {
  ampEvnt_list = append(ampEvnt_list, ae)
  aeCount_list = append(aeCount_list, sum(clin_amplification_therapeutics_ggplot_rmNone$ampEvnt==ae))
}

ampEvnt_vec <- data.frame(matrix(unlist(ampEvnt_list),byrow=T),stringsAsFactors=FALSE)
aeCount_vec <- data.frame(matrix(unlist(aeCount_list),byrow=T),stringsAsFactors=FALSE)

temp = cbind(ampEvnt_vec, aeCount_vec)
colnames(temp)[1] = "ampEvnt"
colnames(temp)[2] = "aeCount"
temp = temp[order(temp$aeCount, decreasing=T),]

# RESPONSIVE VARIANTS #
temp = temp[1:11,]

# RESISTANT VARIANTS #
#temp = temp[1:6,]

for (ae in unique(clin_amplification_therapeutics_ggplot_rmNone$ampEvnt)) {
  for (ae2 in unique(temp$ampEvnt)) {
    if (ae == ae2) {
      
      clin_amplification_therapeutics_ggplot_rmNone$tempCol[clin_amplification_therapeutics_ggplot_rmNone$ampEvnt==ae] = "XXXXX"
      
    }
  }
}

clin_amplification_therapeutics_ggplot_rmNone$ampEvnt[is.na(clin_amplification_therapeutics_ggplot_rmNone$tempCol)] = "other"
clin_amplification_therapeutics_ggplot_rmNone$tempCol = NULL

# Responsive variants #
clin_amplification_therapeutics_ggplot_rmNone$ampEvnt = factor(clin_amplification_therapeutics_ggplot_rmNone$ampEvnt, levels=c(
"CCND1:AMP", "ERBB2:AMP", "FGFR1:AMP", "EGFR:AMP", "CCND2:AMP", "CDK4:AMP",
"MDM2:AMP", "KRAS:AMP", "MET:AMP", "CDK6:AMP", "FGFR2:AMP", "other"))

# Resistant variants #
# clin_amplification_therapeutics_ggplot_rmNone$ampEvnt = factor(clin_amplification_therapeutics_ggplot_rmNone$ampEvnt, levels=c(
# "CCND1:AMP", "ERBB2:AMP", "MDM2:AMP", "BRAF:AMP", "KIT:AMP", "MET:AMP", "other"))

library(scales)
library(RColorBrewer)
source("../global_aes_out.R")

################################ RESPONSIVE FIGURES ################################ 

###GGPLOT DRUG PERCENT###
p = ggplot(data=clin_amplification_therapeutics_ggplot, aes(x = plot_age, alpha = ab_drug))
p = p + facet_wrap( .~acronym, drop=T,scale="fixed", nrow=1)
p = p + geom_bar(position = "fill", fill = "dodgerblue4")
p = p + labs(y = "cases with treatment option(s) (%)", alpha = "highest predicted\nactionability level") + theme_bw()
p = p + scale_y_continuous(labels=percent)
p = p + scale_alpha_manual(values=c("A on-label" = 1, "A off-label" = 0.7, "B on-label" = 0.35, "B off-label" = 0.1, "None" = 0),
                           breaks=c("A on-label", "A off-label", "B on-label", "B off-label"))
p = p + theme(legend.position = "top")
p = p + xlab(element_blank())
p = p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p
fn = "out/YoungAdult_Amplification_ClinAct_BarplotBinary_Percent.pdf"
ggsave(fn, w=7, h=4, useDingbats=FALSE) ###

###GGPLOT DRUG AMPEVNT PERCENT###
p = ggplot(data=clin_amplification_therapeutics_ggplot_rmNone, aes(x = plot_age, fill = ampEvnt))
p = p + facet_wrap( .~acronym, drop=T,scale="fixed", nrow=1)
p = p + geom_bar(position = "fill")
p = p + scale_fill_brewer(palette = "Paired")
p = p + labs(y = "clinically druggable amp. events (%)", x = "age at pathologic diagnosis (yrs.)",
             fill = "events with A\nor B level drug(s)") + theme_bw() ###
p = p + scale_y_reverse(labels = scales::percent)
p = p + theme(strip.background = element_blank(),strip.text.x = element_blank())
p = p + theme(legend.position = "bottom")
p
fn = "out/YoungAdult_AmplificationDrug_AmpEvnt_BarplotBinary_Percent.pdf"
ggsave(fn, w=7, h=4.5, useDingbats=FALSE)

###GGPLOT DRUG COUNT###
p = ggplot(data=clin_amplification_therapeutics_ggplot, aes(x = plot_age, fill = ab_drug))
p = p + facet_wrap( .~acronym, drop=T,scale="fixed", nrow=1)
p = p + geom_bar(position = "stack") 
p = p + labs(fill = "highest predicted\nactionability level", y = "individual cases",
             x = "age at pathologic diagnosis (yrs.)") + theme_bw()
p = p + scale_fill_manual(values=c("None" = "lightslategray", "A on-label" = "dodgerblue1", "A off-label" = "seagreen3",
                                   "B on-label" = "lightgoldenrod", "B off-label" = "lightcoral"),
                          breaks=c("A on-label", "A off-label", "B on-label", "B off-label", "None"))
p = p + theme(legend.position="top")
p
fn = "out/YoungAdult_Amplification_ClinAct_BarplotBinary_Count.pdf"
ggsave(fn, w=7, h=3.5, useDingbats=FALSE)

################################ RESISTANT FIGURES ################################ 

###GGPLOT RESISTANCE PERCENT###
p = ggplot(data=clin_amplification_therapeutics_ggplot, aes(x = plot_age, alpha = ab_drug))
p = p + facet_wrap( .~acronym, drop=T,scale="fixed", nrow=1)
p = p + geom_bar(position = "fill", fill = "red")
p = p + labs(y = "cases resistant to treatment(s) (%)", alpha = "highest\npredicted\nresistance\nlevel") + theme_bw()
p = p + scale_y_continuous(labels=percent)
p = p + scale_alpha_manual(values=c("A on-label" = 1, "A off-label" = 0.7, "B on-label" = 0.35, "B off-label" = 0.1, "None" = 0),
                           breaks=c("A on-label", "A off-label", "B on-label", "B off-label"))
p = p + theme(legend.position = "right")
p = p + xlab(element_blank())
p = p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p
fn = "out/YoungAdult_Amplification_ClinResist_BarplotBinary_Percent.pdf"
ggsave(fn, w=5.25, h=3.5, useDingbats=FALSE) ###

###GGPLOT RESISTANCE GENEVAR PERCENT###
p = ggplot(data=clin_amplification_therapeutics_ggplot_rmNone, aes(x = plot_age, fill = ampEvnt))
p = p + facet_wrap( .~acronym, drop=T,scale="fixed", nrow=1)
p = p + geom_bar(position = "fill")
p = p + scale_fill_brewer(palette = "Paired")
p = p + labs(y = "clinically resistant gene variants (%)", x = "age at pathologic diagnosis (yrs.)",
             fill = "variants\nresistant\nto A or B\nlevel drug(s)") + theme_bw() ###
p = p + scale_y_reverse(labels = scales::percent)
p = p + theme(strip.background = element_blank(),strip.text.x = element_blank())
p = p + theme(legend.position = "right")
p
fn = "out/YoungAdult_AmplificationResist_AmpEvnt_BarplotBinary_Percent.pdf"
ggsave(fn, w=5.25, h=3.5, useDingbats=FALSE)

###GGPLOT RESISTANCE COUNT###
p = ggplot(data=clin_amplification_therapeutics_ggplot, aes(x = plot_age, fill = ab_drug))
p = p + facet_wrap( .~acronym, drop=T,scale="fixed", nrow=1)
p = p + geom_bar(position = "stack") 
p = p + labs(fill = "highest predicted\nresistance level", y = "individual cases",
             x = "age at pathologic diagnosis (yrs.)") + theme_bw()
p = p + scale_fill_manual(values=c("None" = "lightslategray", "A on-label" = "dodgerblue1", "A off-label" = "seagreen3",
                                   "B on-label" = "lightgoldenrod", "B off-label" = "lightcoral"),
                          breaks=c("A on-label", "A off-label", "B on-label", "B off-label", "None"))
p = p + theme(legend.position="top")
p
fn = "out/YoungAdult_Amplification_ClinResist_BarplotBinary_Count.pdf"
ggsave(fn, w=4.5, h=3.5, useDingbats=FALSE)

################################ RESPONSIVE SUMMARY TABLES ################################

library(kableExtra)

ya_amplification_summ_df <- data.frame(matrix(ncol = 7, nrow = 0))
lo_amplification_summ_df <- data.frame(matrix(ncol = 7, nrow = 0))

for (cancer in unique(clin_amplification_therapeutics_ggplot$acronym)) {
  
  temp = clin_amplification_therapeutics_ggplot[clin_amplification_therapeutics_ggplot$acronym==cancer,]
  
  ya_sample_size = sum(temp$age_binary==TRUE)
  lo_sample_size = sum(temp$age_binary==FALSE)
  
  ya_AonLab = sum(temp$ab_drug=="A on-label" & temp$age_binary==TRUE)
  ya_AonLab_perc = round((ya_AonLab/ya_sample_size)*100, digits=1)
  
  ya_AoffLab = sum(temp$ab_drug=="A off-label" & temp$age_binary==TRUE)
  ya_AoffLab_perc = round((ya_AoffLab/ya_sample_size)*100, digits=1)
  
  ya_BonLab = sum(temp$ab_drug=="B on-label" & temp$age_binary==TRUE)
  ya_BonLab_perc = round((ya_BonLab/ya_sample_size)*100, digits=1)
  
  ya_BoffLab = sum(temp$ab_drug=="B off-label" & temp$age_binary==TRUE)
  ya_BoffLab_perc = round((ya_BoffLab/ya_sample_size)*100, digits=1)
  
  ya_None = sum(temp$ab_drug=="None" & temp$age_binary==TRUE)
  ya_None_perc = round((ya_None/ya_sample_size)*100, digits=1)
  
  lo_AonLab = sum(temp$ab_drug=="A on-label" & temp$age_binary==FALSE)
  lo_AonLab_perc = round((lo_AonLab/lo_sample_size)*100, digits=1)
  
  lo_AoffLab = sum(temp$ab_drug=="A off-label" & temp$age_binary==FALSE)
  lo_AoffLab_perc = round((lo_AoffLab/lo_sample_size)*100, digits=1)
  
  lo_BonLab = sum(temp$ab_drug=="B on-label" & temp$age_binary==FALSE)
  lo_BonLab_perc = round((lo_BonLab/lo_sample_size)*100, digits=1)
  
  lo_BoffLab = sum(temp$ab_drug=="B off-label" & temp$age_binary==FALSE)
  lo_BoffLab_perc = round((lo_BoffLab/lo_sample_size)*100, digits=1)
  
  lo_None = sum(temp$ab_drug=="None" & temp$age_binary==FALSE)
  lo_None_perc = round((lo_None/lo_sample_size)*100, digits=1)
  
  ya_row = cbind(cancer, ya_sample_size, ya_AonLab_perc, ya_AoffLab_perc, ya_BonLab_perc, ya_BoffLab_perc, ya_None_perc)
  lo_row = cbind(cancer, lo_sample_size, lo_AonLab_perc, lo_AoffLab_perc, lo_BonLab_perc, lo_BoffLab_perc, lo_None_perc)
  
  ya_amplification_summ_df <- rbind(ya_amplification_summ_df, ya_row)
  lo_amplification_summ_df <- rbind(lo_amplification_summ_df, lo_row)
  
}

##############################################################################################################################

ya_amplification_summ_df$cancer = as.character(ya_amplification_summ_df$cancer)
ya_amplification_summ_df = ya_amplification_summ_df[order(ya_amplification_summ_df$cancer),]
rownames(ya_amplification_summ_df) <- ya_amplification_summ_df$cancer
ya_amplification_summ_df$cancer = NULL

colnames(ya_amplification_summ_df) <- c("Young adult cases","% A on-label","% A off-label","% B on-label","% B off-label", "% None")

young_adult_table = kable(ya_amplification_summ_df,"html",align="l") %>% kable_styling("striped",full_width = FALSE)

young_adult_table

##############################################################################################################################

lo_amplification_summ_df$cancer = as.character(lo_amplification_summ_df$cancer)
lo_amplification_summ_df = lo_amplification_summ_df[order(lo_amplification_summ_df$cancer),]
rownames(lo_amplification_summ_df) <- lo_amplification_summ_df$cancer
lo_amplification_summ_df$cancer = NULL

colnames(lo_amplification_summ_df) <- c("Later onset cases","% A on-label","% A off-label","% B on-label","% B off-label", "% None")

later_onset_table = kable(lo_amplification_summ_df,"html",align="l") %>% kable_styling("striped",full_width = FALSE)

later_onset_table

################################ RESPONSIVE AMPEVNT TABLES ################################

ampEvnt_list = list(); aeCount_list = list()

for (ae in unique(clin_amplification_actionable$ampEvnt)) {
  ampEvnt_list = append(ampEvnt_list, ae)
  aeCount_list = append(aeCount_list, sum(clin_amplification_actionable$ampEvnt==ae))
}

ampEvnt_vec <- data.frame(matrix(unlist(ampEvnt_list),byrow=T),stringsAsFactors=FALSE)
aeCount_vec <- data.frame(matrix(unlist(aeCount_list),byrow=T),stringsAsFactors=FALSE)

temp = cbind(ampEvnt_vec, aeCount_vec)
colnames(temp)[1] = "ampEvnt"
colnames(temp)[2] = "aeCount"
temp = temp[order(temp$aeCount, decreasing=T),]

temp = temp[1:10,]

###############################################

ya_summ_ampEvnt_df <- data.frame(matrix(ncol = 12, nrow = 0))
lo_summ_ampEvnt_df <- data.frame(matrix(ncol = 12, nrow = 0))

young_adult_n = data.frame(matrix(ncol = 1, nrow = 0))
later_onset_n = data.frame(matrix(ncol = 1, nrow = 0))

amp_evnt_vec = c("CCND1:AMP", "FGFR1:AMP", "ERBB2:AMP", "EGFR:AMP", "MDM2:AMP",
                 "CCND2:AMP", "KRAS:AMP", "CDK4:AMP", "CDK6:AMP", "MET:AMP")

for (cancer in unique(clin_amplification_therapeutics_all$acronym)) {

  current_cancer_df = clin_amplification_therapeutics_all[clin_amplification_therapeutics_all$acronym==cancer,]
  
  unique_cases = current_cancer_df[!duplicated(current_cancer_df$bcr_patient_barcode),]
  ya_sample_size = sum(unique_cases$age_binary==TRUE)
  lo_sample_size = sum(unique_cases$age_binary==FALSE)
  
  ya_ampEvnt_row <- data.frame(matrix(ncol = 0, nrow = 1))
  lo_ampEvnt_row <- data.frame(matrix(ncol = 0, nrow = 1))
  
  for (ampEvnt in amp_evnt_vec) {
    
    num_current_ampEvnt_yng = sum(current_cancer_df$age_binary==TRUE & current_cancer_df$ampEvnt==ampEvnt)
    num_current_ampEvnt_ltr = sum(current_cancer_df$age_binary==FALSE & current_cancer_df$ampEvnt==ampEvnt)
    
    ya_current_ampEvnt_perc = round((num_current_ampEvnt_yng/ya_sample_size)*100, digits=1)
    lo_current_ampEvnt_perc = round((num_current_ampEvnt_ltr/lo_sample_size)*100, digits=1)
    
    ya_ampEvnt_row = cbind(ya_ampEvnt_row, ya_current_ampEvnt_perc)
    lo_ampEvnt_row = cbind(lo_ampEvnt_row, lo_current_ampEvnt_perc)
    
  }
  
  young_adult_n = rbind(young_adult_n, ya_sample_size)
  later_onset_n = rbind(later_onset_n, lo_sample_size)
  
  ya_ampEvnt_row = cbind(ya_ampEvnt_row, cancer)
  lo_ampEvnt_row = cbind(lo_ampEvnt_row, cancer)
  
  ya_summ_ampEvnt_df <- rbind(ya_summ_ampEvnt_df, ya_ampEvnt_row)
  lo_summ_ampEvnt_df <- rbind(lo_summ_ampEvnt_df, lo_ampEvnt_row)
  
}

##############################################################################################################################

colnames(ya_summ_ampEvnt_df) <- c("% CCND1:AMP", "% FGFR1:AMP", "% ERBB2:AMP", "% EGFR:AMP", "% MDM2:AMP",
                                  "% CCND2:AMP", "% KRAS:AMP", "% CDK4:AMP", "% CDK6:AMP", "% MET:AMP", "cancer")

ya_summ_ampEvnt_df = add_column(ya_summ_ampEvnt_df, young_adult_n, .before ="% CCND1:AMP")

ya_summ_ampEvnt_df$cancer = as.character(ya_summ_ampEvnt_df$cancer)
ya_summ_ampEvnt_df = ya_summ_ampEvnt_df[order(ya_summ_ampEvnt_df$cancer),]
rownames(ya_summ_ampEvnt_df) <- ya_summ_ampEvnt_df$cancer
ya_summ_ampEvnt_df$cancer = NULL

colnames(ya_summ_ampEvnt_df)[1] = "Young adult cases"

young_ampEvnt_table = kable(ya_summ_ampEvnt_df,"html",align="l") %>% kable_styling("striped",full_width = FALSE)

young_ampEvnt_table

##############################################################################################################################

colnames(lo_summ_ampEvnt_df) <- c("% CCND1:AMP", "% FGFR1:AMP", "% ERBB2:AMP", "% EGFR:AMP", "% MDM2:AMP",
                                  "% CCND2:AMP", "% KRAS:AMP", "% CDK4:AMP", "% CDK6:AMP", "% MET:AMP", "cancer")

lo_summ_ampEvnt_df = add_column(lo_summ_ampEvnt_df, later_onset_n, .before ="% CCND1:AMP")

lo_summ_ampEvnt_df$cancer = as.character(lo_summ_ampEvnt_df$cancer)
lo_summ_ampEvnt_df = lo_summ_ampEvnt_df[order(lo_summ_ampEvnt_df$cancer),]
rownames(lo_summ_ampEvnt_df) <- lo_summ_ampEvnt_df$cancer
lo_summ_ampEvnt_df$cancer = NULL

colnames(lo_summ_ampEvnt_df)[1] = "Later onset cases"

later_ampEvnt_table = kable(lo_summ_ampEvnt_df,"html",align="l") %>% kable_styling("striped",full_width = FALSE)

later_ampEvnt_table
