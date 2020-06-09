##### YoungAdult_Somatic_Clin_Assoc.R #####
# William Lee @ January 2020
# Updated June 2020

### analytic pipeline ###

library(tidyverse)

# set working directory for dev and debug
bdir = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/clinical_associations"
setwd(bdir)

# civic_raw_data updated 06/04/2020
civic_raw_data = "01-May-2020-ClinicalEvidenceSummaries.tsv" # 3431 obs. of 41 variables
# civic_raw_data = "01-Apr-2020-ClinicalEvidenceSummaries.tsv" # 3393 obs. of 41 variables
civic_df = read.table(header=T, quote = "", sep="\t", fill=T, file = civic_raw_data, stringsAsFactors=FALSE)

### civic data preparation 
civic_df$geneVar = paste(civic_df$gene,civic_df$variant,sep=":")
civic_df$doid[is.na(civic_df$doid)] = "XXXXX"
### end of civic data preparation 

# cgi_raw_data updated 06/04/2020
cgi_raw_data = "cgi_biomarkers_per_variant.tsv" # 1442 obs. of 28 variables
cgi_df = read.table(header=T, quote = "", sep="\t", fill=T, file = cgi_raw_data, stringsAsFactors=FALSE)

### cgi data preparation 
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
### end of cgi data preparation 

# onco_raw_data updated 06/04/2020
onco_raw_data = "oncokb_biomarker_drug_associations.tsv"
onco_df = read.table(header=T, quote = "", sep="\t", fill=T, file = onco_raw_data, stringsAsFactors=FALSE)

### onco kb data preparation
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

onco_df$Alterations = str_replace_all(onco_df$Alterations, fixed(" "), "")

onco_df_fin <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(onco_df_fin) <- colnames(onco_df)
i = 1
for (element in onco_df$Alterations) {
  for (mut in (unlist(strsplit(element, ",")))) {
    curRow = onco_df[i,]; curRow$Alterations = mut
    onco_df_fin = rbind(onco_df_fin, curRow)
  }
  i = i + 1
}

onco_df_fin$geneVar = paste(onco_df_fin$Gene,onco_df_fin$Alterations,sep=":")

onco_df_fin$clin.sig[onco_df_fin$Level=="1" | onco_df_fin$Level=="2" | onco_df_fin$Level=="3"] = "Responsive"
onco_df_fin$clin.sig[onco_df_fin$Level=="R1" | onco_df_fin$Level=="R2"] = "Resistant"
### end of onco kb data preparation

##### the following code generates combination df fed into analytic pipeline ##### 
doid1 = civic_df$doid
doid2 = cgi_df$doid
doid3 = onco_df_fin$doid
doid = c(doid1, doid2, doid3)

drug1 = civic_df$drugs
drug2 = cgi_df$Drug
drug3 = onco_df_fin$Drugs
drugs = c(drug1, drug2, drug3)

geneVar1 = civic_df$geneVar
geneVar2 = cgi_df$individual_mutation
geneVar3 = onco_df_fin$geneVar
geneVar = c(geneVar1, geneVar2, geneVar3)

elvl1 = civic_df$evidence_level
elvl2 = cgi_df$Evidence.level
elvl3 = onco_df_fin$Evidence.Level
evidence_level = c(elvl1, elvl2, elvl3)

clin1 = civic_df$clinical_significance
clin2 = cgi_df$Association
clin3 = onco_df_fin$clin.sig
clinical_significance = c(clin1, clin2, clin3)

poten_clin_act_3db = as.data.frame(cbind(doid, drugs, geneVar, evidence_level, clinical_significance))

######################################################################################################

### MC3 mutation call file (only include the likely driver for the first-pass analysis) ###
somatic_likelyfunctional_driver_f = "~/Box/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/MC3_Ellrott_CellSys2018/mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample_likelyDriver.tsv"
somatic_likelyfunctional_driver = read.table(header=T, quote = "", sep="\t", file = somatic_likelyfunctional_driver_f, stringsAsFactors=FALSE)

somatic_f = "~/Box/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/MC3_Ellrott_CellSys2018/mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample.gz"
somatic = read.table(header=T, quote = "", sep="\t", file = gzfile(somatic_f), stringsAsFactors=FALSE)
unique_somatic_samples = unique(substr(somatic$Tumor_Sample_Barcode,1,12))

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
clin_merge_subtype_avail = clin_merge_subtype[clin_merge_subtype$bcr_patient_barcode %in% unique_somatic_samples,]
clin_merge_subtype_avail_complete = clin_merge_subtype_avail[complete.cases(clin_merge_subtype_avail[,2:4]) & clin_merge_subtype_avail$gender != "",]

# quick check of merged data
cat("TCGA cases count available for this analysis\n")
nrow(clin_merge_subtype_avail_complete)
cat("Cancer type distribution of TCGA cases with young onset cancer\n")
table(clin_merge_subtype_avail_complete$age_at_initial_pathologic_diagnosis<=50,clin_merge_subtype_avail_complete$acronym)

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

unique_genevar_3db = unique(poten_clin_act_3db$geneVar)

clin_somatic_therapeutics <- data.frame(matrix(ncol = 35, nrow = 0))
colnames(clin_somatic_therapeutics) <- c("bcr_patient_barcode","acronym","age_at_initial_pathologic_diagnosis","gender",
                                         "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","SAMPLE_BARCODE",
                                         "DISEASE","SUBTYPE","age_binary","doid","Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode","HGVSp_Short","geneVar","ab_drug")

# conduct tests by each cancer (denoted as acronym)
for (cancer in unique(clin_merge_subtype_avail_complete$acronym)){
  
  if (cancer == "BRCA" | cancer == "CESC" | cancer == "COAD" | cancer == "HNSC" |
      cancer == "KIRC" | cancer == "KIRP" | cancer == "LGG"  | cancer == "LIHC" |
      cancer == "OV"   | cancer == "PCPG" | cancer == "SARC" | cancer == "SKCM" | 
      cancer == "THCA" | cancer == "UCEC") {
    
    clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete[clin_merge_subtype_avail_complete$acronym==cancer,]
    clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete_c[complete.cases(clin_merge_subtype_avail_complete_c),]
    
    # conduct the test by gene
    for (gene in unique(somatic_likelyfunctional_driver$Hugo_Symbol)){
      
      # subset mutation to that gene and get unique mutation
      somatic_likelyfunctional_driver_g = somatic_likelyfunctional_driver[somatic_likelyfunctional_driver$Hugo_Symbol==gene,]
      somatic_likelyfunctional_driver_g_uniq = somatic_likelyfunctional_driver_g[!duplicated(somatic_likelyfunctional_driver_g$bcr_patient_barcode),]
      clin_merge_subtype_avail_complete_c_merge_somatic = merge(clin_merge_subtype_avail_complete_c,somatic_likelyfunctional_driver_g_uniq, by= "bcr_patient_barcode", all.x=T, all.y=F)
      
      clin_merge_subtype_avail_complete_c_merge_somatic = clin_merge_subtype_avail_complete_c_merge_somatic[complete.cases(clin_merge_subtype_avail_complete_c_merge_somatic$Hugo_Symbol),]
      clin_merge_subtype_avail_complete_c_merge_somatic = clin_merge_subtype_avail_complete_c_merge_somatic[complete.cases(clin_merge_subtype_avail_complete_c_merge_somatic$HGVSp_Short),]
      clin_merge_subtype_avail_complete_c_merge_somatic$HGVSp_Short = str_remove(clin_merge_subtype_avail_complete_c_merge_somatic$HGVSp_Short, "p.")
      clin_merge_subtype_avail_complete_c_merge_somatic$geneVar = paste(clin_merge_subtype_avail_complete_c_merge_somatic$Hugo_Symbol,clin_merge_subtype_avail_complete_c_merge_somatic$HGVSp_Short,sep=":")
      
      clin_merge_no_match = clin_merge_subtype_avail_complete_c_merge_somatic[!(clin_merge_subtype_avail_complete_c_merge_somatic$geneVar %in% unique_genevar_3db),]
      clin_merge_subtype_avail_complete_c_merge_somatic = clin_merge_subtype_avail_complete_c_merge_somatic[clin_merge_subtype_avail_complete_c_merge_somatic$geneVar %in% unique_genevar_3db,]
      
      for (uniq_geneVar in unique(clin_merge_subtype_avail_complete_c_merge_somatic$geneVar)) {
        
        # poten_clin_act_som contains clin_merge_subtype_avail_complete_c_merge_somatic rows that have gene variant of interest 
        poten_clin_act_som = clin_merge_subtype_avail_complete_c_merge_somatic[clin_merge_subtype_avail_complete_c_merge_somatic$geneVar==uniq_geneVar,]
        
        # poten_clin_act_var contains only rows that have gene variant of interest 
        poten_clin_act_var = poten_clin_act_3db[poten_clin_act_3db$geneVar==uniq_geneVar,]
        
        # replaces blank cells in poten_clin_act_var with NA (so complete.cases() will work)
        poten_clin_act_var[poten_clin_act_var == ""] <- NA
        
        # confirm that all rows in poten_clin_act_var have evidence levels and drugs
        poten_clin_act_var = poten_clin_act_var[complete.cases(poten_clin_act_var$evidence_level),]
        poten_clin_act_var = poten_clin_act_var[complete.cases(poten_clin_act_var$drugs),]
        
        # skips to next loop iteration if poten_clin_act_var has no obs.
        if (dim(poten_clin_act_var)[1] == 0) {
          poten_clin_act_som$ab_drug = ""
          clin_somatic_therapeutics <- rbind(clin_somatic_therapeutics, poten_clin_act_som)
          next
        }
        
        # instantiates empty list (elements to be appended)
        current_list = list()
        
        # initiates counter at 1
        i = 1
        
        for (elvl in poten_clin_act_var$evidence_level) {
          
          if (elvl == "A" & (poten_clin_act_var[i, "clinical_significance"] == "Sensitivity/Response" |
                             poten_clin_act_var[i, "clinical_significance"] == "Responsive")) {
            
            # if (elvl == "A" & (poten_clin_act_var[i, "clinical_significance"] == "Resistance" |
            #                    poten_clin_act_var[i, "clinical_significance"] == "Resistant")) {
            
            if (poten_clin_act_var[i, "doid"] == unique(poten_clin_act_som$doid)) {
              current_drug = paste(poten_clin_act_var[i, "drugs"], "A on-label", sep=": ")
              current_list = append(current_list, current_drug)
            }
            
            else {
              current_drug = paste(poten_clin_act_var[i, "drugs"], "A off-label", sep=": ")
              current_list = append(current_list, current_drug)
            }
            
            i = i + 1
            
          }
          
          else if (elvl == "B" & (poten_clin_act_var[i, "clinical_significance"] == "Sensitivity/Response" |
                                  poten_clin_act_var[i, "clinical_significance"] == "Responsive")) {
            
            # else if (elvl == "B" & (poten_clin_act_var[i, "clinical_significance"] == "Resistance" |
            #                         poten_clin_act_var[i, "clinical_significance"] == "Resistant")) {
            
            if (poten_clin_act_var[i, "doid"] == unique(poten_clin_act_som$doid)) {
              current_drug = paste(poten_clin_act_var[i, "drugs"], "B on-label", sep=": ")
              current_list = append(current_list, current_drug)
            }
            
            else {
              current_drug = paste(poten_clin_act_var[i, "drugs"], "B off-label", sep=": ")
              current_list = append(current_list, current_drug)
            }
            
            i = i + 1
            
          }
          
          else {i = i + 1}
          
        }
        
        current_list = unique(current_list)
        therapeutic_drugs = paste(current_list, collapse = " & ")
        poten_clin_act_som$ab_drug = therapeutic_drugs
        
        clin_somatic_therapeutics <- rbind(clin_somatic_therapeutics, poten_clin_act_som)
        
        
      }
      
      if (dim(clin_merge_no_match)[1] != 0) {
        clin_merge_no_match$ab_drug = ""
        clin_somatic_therapeutics <- rbind(clin_somatic_therapeutics, clin_merge_no_match)
      }
    }
  }
}

# RESPONSIVE VARIANTS #
# clin_somatic_therapeutics: 17523 obs. before removing duplicates

# RESISTANT VARIANTS #
# clin_somatic_therapeutics: 17523 obs. before removing duplicates

duplicated_patient_samples = clin_somatic_therapeutics[duplicated(clin_somatic_therapeutics$bcr_patient_barcode),1]

for (duplicate in duplicated_patient_samples) {
  
  current_duplicates_df = clin_somatic_therapeutics[clin_somatic_therapeutics$bcr_patient_barcode==duplicate,]
  clin_somatic_therapeutics = clin_somatic_therapeutics[clin_somatic_therapeutics$bcr_patient_barcode!=duplicate,]
  
  drugs_as_string = paste(current_duplicates_df$ab_drug, collapse = " ")
  
  if (str_detect(drugs_as_string, "A on-label")) {
    current_duplicates_df = current_duplicates_df[str_detect(current_duplicates_df$ab_drug, "A on-label"),]
    clin_somatic_therapeutics = rbind(clin_somatic_therapeutics, current_duplicates_df[1,])
  }
  
  else if (str_detect(drugs_as_string, "A off-label")) {
    current_duplicates_df = current_duplicates_df[str_detect(current_duplicates_df$ab_drug, "A off-label"),]
    clin_somatic_therapeutics = rbind(clin_somatic_therapeutics, current_duplicates_df[1,])
  }
  
  else if (str_detect(drugs_as_string, "B on-label")) {
    current_duplicates_df = current_duplicates_df[str_detect(current_duplicates_df$ab_drug, "B on-label"),]
    clin_somatic_therapeutics = rbind(clin_somatic_therapeutics, current_duplicates_df[1,])
  }
  
  else if (str_detect(drugs_as_string, "B off-label")) {
    current_duplicates_df = current_duplicates_df[str_detect(current_duplicates_df$ab_drug, "B off-label"),]
    clin_somatic_therapeutics = rbind(clin_somatic_therapeutics, current_duplicates_df[1,])
  }
  
  else {
    clin_somatic_therapeutics = rbind(clin_somatic_therapeutics, current_duplicates_df[1,])
  }
  
}

# RESPONSIVE VARIANTS #
# clin_somatic_therapeutics: 4665 obs. after removing duplicates

# RESISTANT VARIANTS #
# clin_somatic_therapeutics: 4665 obs. after removing duplicates

### the following code produces the ggplot-ready dataframe ###

# replaces blank cells in clin_somatic_therapeutics with NA   
clin_somatic_therapeutics[clin_somatic_therapeutics == ""] <- NA

clin_somatic_therapeutics_ggplot <- data.frame(matrix(ncol = 35, nrow = 0))
colnames(clin_somatic_therapeutics_ggplot) <- colnames(clin_somatic_therapeutics)

for (cancer in unique(clin_somatic_therapeutics$acronym)){
  
  i = 1
  
  clin_somatic_therapeutics_c = clin_somatic_therapeutics[clin_somatic_therapeutics$acronym==cancer,]
  
  total_reg_cases = sum(clin_somatic_therapeutics_c$age_binary==FALSE & !is.na(clin_somatic_therapeutics_c$ab_drug))
  total_yng_cases = sum(clin_somatic_therapeutics_c$age_binary==TRUE & !is.na(clin_somatic_therapeutics_c$ab_drug))
  
  if (total_reg_cases >= 10 & total_yng_cases >= 10) {
    
    for (element in clin_somatic_therapeutics_c$ab_drug) {
      
      if (is.na(element)) {i = i + 1}
      
      else if (str_detect(element, "A on-label")) {clin_somatic_therapeutics_c[i, "ab_drug"] =  "A on-label"; i = i + 1}
      
      else if (str_detect(element, "A off-label")) {clin_somatic_therapeutics_c[i, "ab_drug"] =  "A off-label"; i = i + 1}
      
      else if (str_detect(element, "B on-label")) {clin_somatic_therapeutics_c[i, "ab_drug"] =  "B on-label"; i = i + 1}
      
      else if (str_detect(element, "B off-label")) {clin_somatic_therapeutics_c[i, "ab_drug"] =  "B off-label"; i = i + 1}
    }
    
    clin_somatic_therapeutics_ggplot <- rbind(clin_somatic_therapeutics_ggplot, clin_somatic_therapeutics_c)
    
  }
}

# RESPONSIVE VARIANTS #
# clin_somatic_therapeutics_ggplot: 3062 obs. 

# RESISTANT VARIANTS #
# clin_somatic_therapeutics_ggplot: ? obs. 

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

clin_somatic_therapeutics_ggplot$ab_drug[is.na(clin_somatic_therapeutics_ggplot$ab_drug)] = "None"

clin_somatic_therapeutics_ggplot_rmNone = clin_somatic_therapeutics_ggplot[clin_somatic_therapeutics_ggplot$ab_drug!="None",]
clin_somatic_therapeutics_ggplot_None = clin_somatic_therapeutics_ggplot[clin_somatic_therapeutics_ggplot$ab_drug=="None",]

clin_somatic_therapeutics_ggplot$ab_drug = factor(clin_somatic_therapeutics_ggplot$ab_drug, levels=c("None", "B off-label", "B on-label", "A off-label", "A on-label"))

clin_somatic_therapeutics_ggplot$plot_age[clin_somatic_therapeutics_ggplot$age_binary==TRUE] = "<= 50"
clin_somatic_therapeutics_ggplot$plot_age[clin_somatic_therapeutics_ggplot$age_binary==FALSE] = "> 50"

clin_somatic_therapeutics_ggplot_rmNone$plot_age[clin_somatic_therapeutics_ggplot_rmNone$age_binary==TRUE] = "<= 50"
clin_somatic_therapeutics_ggplot_rmNone$plot_age[clin_somatic_therapeutics_ggplot_rmNone$age_binary==FALSE] = "> 50"

clin_somatic_therapeutics_ggplot_None$plot_age[clin_somatic_therapeutics_ggplot_None$age_binary==TRUE] = "<= 50"
clin_somatic_therapeutics_ggplot_None$plot_age[clin_somatic_therapeutics_ggplot_None$age_binary==FALSE] = "> 50"

### Selects for top gene variants ###

geneVar_list = list(); gvCount_list = list()

for (gv in unique(clin_somatic_therapeutics_ggplot_rmNone$geneVar)) {
  geneVar_list = append(geneVar_list, gv)
  gvCount_list = append(gvCount_list, sum(clin_somatic_therapeutics_ggplot_rmNone$geneVar==gv))
}

geneVar_vec <- data.frame(matrix(unlist(geneVar_list),byrow=T),stringsAsFactors=FALSE)
gvCount_vec <- data.frame(matrix(unlist(gvCount_list),byrow=T),stringsAsFactors=FALSE)

temp = cbind(geneVar_vec, gvCount_vec)
colnames(temp)[1] = "geneVar"
colnames(temp)[2] = "gvCount"
temp = temp[order(temp$gvCount, decreasing=T),]
temp = temp[1:11,]

for (gv in unique(clin_somatic_therapeutics_ggplot_rmNone$geneVar)) {
  for (gv2 in unique(temp$geneVar)) {
    if (gv == gv2) {
      
      clin_somatic_therapeutics_ggplot_rmNone$tempCol[clin_somatic_therapeutics_ggplot_rmNone$geneVar==gv] = "XXXXX"
      
    }
  }
}

clin_somatic_therapeutics_ggplot_rmNone$geneVar[is.na(clin_somatic_therapeutics_ggplot_rmNone$tempCol)] = "other"
clin_somatic_therapeutics_ggplot_rmNone$tempCol = NULL

# Responsive variants #
clin_somatic_therapeutics_ggplot_rmNone$geneVar = factor(clin_somatic_therapeutics_ggplot_rmNone$geneVar, levels=c(
  "BRAF:V600E", "PIK3CA:H1047R", "PIK3CA:E545K", "PIK3CA:E542K", "KRAS:G12V", "KRAS:G13D",
  "AKT1:E17K", "KRAS:G12C", "CTNNB1:S37C", "ERBB3:V104M", "IDH1:R132C", "other"))

# Resistant variants #
# clin_somatic_therapeutics_ggplot_rmNone$geneVar = factor(clin_somatic_therapeutics_ggplot_rmNone$geneVar, levels=c(
# "BRAF:V600E", "PIK3CA:H1047R", "PIK3CA:E545K", "PIK3CA:E542K", "KRAS:G12D", "KRAS:G13D",
# "KRAS:A146T", "KRAS:G12A", "KRAS:G12C", "KRAS:G12S", "NRAS:G12D", "other"))


library(scales)
library(RColorBrewer)
source("../global_aes_out.R")

#
#
#
#
#

################################ RESPONSIVE FIGURES ################################ 

###GGPLOT DRUG PERCENT###
p = ggplot(data=clin_somatic_therapeutics_ggplot, aes(x = plot_age, alpha = ab_drug))
p = p + facet_wrap( .~acronym, drop=T,scale="fixed", nrow=1)
#p = p + facet_grid( .~acronym, drop=T,scale="free")
p = p + geom_bar(position = "fill", fill = "dodgerblue4")
p = p + labs(y = "cases with treatment option(s) (%)", #x = "age at pathologic diagnosis (yrs.)",
             alpha = "highest predicted\nactionability level") + theme_bw()
p = p + scale_y_continuous(labels=percent)
p = p + scale_alpha_manual(values=c("A on-label" = 1, "A off-label" = 0.7, "B on-label" = 0.35, "B off-label" = 0.1, "None" = 0),
                           breaks=c("A on-label", "A off-label", "B on-label", "B off-label"))
p = p + theme(legend.position = "top")
p = p + xlab(element_blank())
p = p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p
fn = "out/YoungAdult_Somatic_ClinAct_BarplotBinary_Percent.pdf"
ggsave(fn, w=7.5, h=4, useDingbats=FALSE) ###

###GGPLOT DRUG GENEVAR PERCENT###
p = ggplot(data=clin_somatic_therapeutics_ggplot_rmNone, aes(x = plot_age, fill = geneVar))
p = p + facet_wrap( .~acronym, drop=T,scale="fixed", nrow=1)
p = p + geom_bar(position = "fill")
p = p + scale_fill_brewer(palette = "Paired")
p = p + labs(y = "clinically druggable gene variants (%)", x = "age at pathologic diagnosis (yrs.)",
             fill = "variants with A\nor B level drug(s)") + theme_bw() ###
p = p + scale_y_reverse(labels = scales::percent)
#p = p + scale_y_continuous(labels=percent)
#p = p + xlab(element_blank()) ###
p = p + theme(strip.background = element_blank(),strip.text.x = element_blank())
p = p + theme(legend.position = "bottom")
p
fn = "out/YoungAdult_SomaticDrug_GeneVar_BarplotBinary_Percent.pdf"
ggsave(fn, w=7.5, h=4.5, useDingbats=FALSE)

###GGPLOT DRUG COUNT###
p = ggplot(data=clin_somatic_therapeutics_ggplot, aes(x = plot_age, fill = ab_drug))
p = p + facet_wrap( .~acronym, drop=T,scale="fixed", nrow=1)
p = p + geom_bar(position = "stack") 
p = p + labs(fill = "highest predicted\nactionability level", y = "individual cases",
             x = "age at pathologic diagnosis (yrs.)") + theme_bw()
p = p + scale_fill_manual(values=c("None" = "lightslategray", "A on-label" = "dodgerblue1", "A off-label" = "seagreen3",
                                   "B on-label" = "lightgoldenrod", "B off-label" = "lightcoral"),
                          breaks=c("A on-label", "A off-label", "B on-label", "B off-label", "None"))
p = p + theme(legend.position="top")
p
fn = "out/YoungAdult_Somatic_ClinAct_BarplotBinary_Count.pdf"
ggsave(fn, w=7.5, h=3.5, useDingbats=FALSE)

################################ RESISTANT FIGURES ################################ 

###GGPLOT RESISTANCE PERCENT###
p = ggplot(data=clin_somatic_therapeutics_ggplot, aes(x = plot_age, alpha = ab_drug))
p = p + facet_wrap( .~acronym, drop=T,scale="fixed", nrow=1)
#p = p + facet_grid( .~acronym, drop=T,scale="free")
p = p + geom_bar(position = "fill", fill = "red")
p = p + labs(y = "cases resistant to treatment(s) (%)", #x = "age at pathologic diagnosis (yrs.)",
             alpha = "highest predicted\nactionability level") + theme_bw()
p = p + scale_y_continuous(labels=percent)
p = p + scale_alpha_manual(values=c("A on-label" = 1, "A off-label" = 0.7, "B on-label" = 0.35, "B off-label" = 0.1, "None" = 0),
                           breaks=c("A on-label", "A off-label", "B on-label", "B off-label"))
p = p + theme(legend.position = "top")
p = p + xlab(element_blank())
p = p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p
fn = "out/YoungAdult_Somatic_ClinResist_BarplotBinary_Percent.pdf"
ggsave(fn, w=7.5, h=4, useDingbats=FALSE) ###

###GGPLOT RESISTANCE GENEVAR PERCENT###
p = ggplot(data=clin_somatic_therapeutics_ggplot_rmNone, aes(x = plot_age, fill = geneVar))
p = p + facet_wrap( .~acronym, drop=T,scale="fixed", nrow=1)
p = p + geom_bar(position = "fill")
p = p + scale_fill_brewer(palette = "Paired")
p = p + labs(y = "clinically resistant gene variants (%)", x = "age at pathologic diagnosis (yrs.)",
             fill = "variants resistant\nto A or B\nlevel drug(s)") + theme_bw() ###
p = p + scale_y_reverse(labels = scales::percent)
#p = p + scale_y_continuous(labels=percent)
#p = p + xlab(element_blank()) ###
p = p + theme(strip.background = element_blank(),strip.text.x = element_blank())
p = p + theme(legend.position = "bottom")
p
fn = "out/YoungAdult_SomaticResist_GeneVar_BarplotBinary_Percent.pdf"
ggsave(fn, w=7.5, h=4.5, useDingbats=FALSE)

###GGPLOT RESISTANCE COUNT###
p = ggplot(data=clin_somatic_therapeutics_ggplot, aes(x = plot_age, fill = ab_drug))
p = p + facet_wrap( .~acronym, drop=T,scale="fixed", nrow=1)
p = p + geom_bar(position = "stack") 
p = p + labs(fill = "highest predicted\nactionability level", y = "individual cases",
             x = "age at pathologic diagnosis (yrs.)") + theme_bw()
p = p + scale_fill_manual(values=c("None" = "lightslategray", "A on-label" = "dodgerblue1", "A off-label" = "seagreen3",
                                   "B on-label" = "lightgoldenrod", "B off-label" = "lightcoral"),
                          breaks=c("A on-label", "A off-label", "B on-label", "B off-label", "None"))
p = p + theme(legend.position="top")
p
fn = "out/YoungAdult_Somatic_ClinResist_BarplotBinary_Count.pdf"
ggsave(fn, w=7.5, h=3.5, useDingbats=FALSE)

################################ FISHER'S EXACT DRUG ANALYSIS ################################ 

fisher_abdrug_df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(fisher_abdrug_df) <- c("cancer", "ab_drug", "p_value")

for (cancer in unique(clin_somatic_therapeutics_ggplot$acronym)) {
  
  current_cancer_df = clin_somatic_therapeutics_ggplot[clin_somatic_therapeutics_ggplot$acronym==cancer,]
  
  for (abdrug in unique(current_cancer_df$ab_drug)) {
    
    num_current_abdrug_ltr = sum(current_cancer_df$age_binary==FALSE & current_cancer_df$ab_drug==abdrug)
    num_current_abdrug_yng = sum(current_cancer_df$age_binary==TRUE & current_cancer_df$ab_drug==abdrug)
    
    num_not_current_abdrug_ltr = sum(current_cancer_df$age_binary==FALSE & current_cancer_df$ab_drug!=abdrug)
    num_not_current_abdrug_yng = sum(current_cancer_df$age_binary==TRUE & current_cancer_df$ab_drug!=abdrug)
    
    yng_set = c(num_not_current_abdrug_yng, num_current_abdrug_yng)
    late_set = c(num_not_current_abdrug_ltr, num_current_abdrug_ltr)
    
    abdrug_fisher = cbind(yng_set, late_set)
    
    abdrug_fisher_pval = fisher.test(abdrug_fisher)$p.val
    
    temp = cbind(cancer, abdrug, abdrug_fisher_pval)
    colnames(temp) <- c("cancer", "ab_drug", "p_value")
    
    fisher_abdrug_df = rbind(fisher_abdrug_df, temp)
    
  }
}

fisher_abdrug_df$p_value = as.numeric(as.character(fisher_abdrug_df$p_value))
fisher_abdrug_df$FDR = p.adjust(fisher_abdrug_df[,"p_value"], method="fdr") 
fisher_abdrug_df=fisher_abdrug_df[order(fisher_abdrug_df$p_value, decreasing=FALSE),]

################################ FISHER'S EXACT GENEVAR ANALYSIS ################################

fisher_genevar_df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(fisher_genevar_df) <- c("cancer", "genevar", "p_value")

for (cancer in unique(clin_somatic_therapeutics_ggplot_rmNone$acronym)) {
  
  current_cancer_df = clin_somatic_therapeutics_ggplot_rmNone[clin_somatic_therapeutics_ggplot_rmNone$acronym==cancer,]
  
  for (genevar in unique(current_cancer_df$geneVar)) {
    
    num_current_genevar_ltr = sum(current_cancer_df$age_binary==FALSE & current_cancer_df$geneVar==genevar)
    num_current_genevar_yng = sum(current_cancer_df$age_binary==TRUE & current_cancer_df$geneVar==genevar)
    
    num_not_current_genevar_ltr = sum(current_cancer_df$age_binary==FALSE & current_cancer_df$geneVar!=genevar)
    num_not_current_genevar_yng = sum(current_cancer_df$age_binary==TRUE & current_cancer_df$geneVar!=genevar)
    
    yng_set = c(num_not_current_genevar_yng, num_current_genevar_yng)
    late_set = c(num_not_current_genevar_ltr, num_current_genevar_ltr)
    
    genevar_fisher = cbind(yng_set, late_set)
    
    genevar_fisher_pval = fisher.test(genevar_fisher)$p.val
    
    temp = cbind(cancer, genevar, genevar_fisher_pval)
    colnames(temp) <- c("cancer", "geneVar", "p_value")
    
    fisher_genevar_df = rbind(fisher_genevar_df, temp)
    
  }
}

fisher_genevar_df$p_value = as.numeric(as.character(fisher_genevar_df$p_value))
fisher_genevar_df$FDR = p.adjust(fisher_genevar_df[,"p_value"], method="fdr") 
fisher_genevar_df=fisher_genevar_df[order(fisher_genevar_df$p_value, decreasing=FALSE),]
