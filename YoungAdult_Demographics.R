##### Demographics_YoungAdult_Assoc.R #####
# William Lee @ March 2020
# work in progress

install.packages("data.table")
install.packages("dplyr")
install.packages("formattable")
install.packages("tidyr")
install.packages("kableExtra")

library(data.table)
library(dplyr)
library(formattable)
library(tidyr)
library(kableExtra)
library(tibble)

# set working directory for dev and debug
bdir = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/Demographics"
setwd(bdir)

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

clin_merge_subtype_avail_complete$cancer = NULL
clin_merge_subtype_avail_complete$ethnicity = NULL
clin_merge_subtype_avail_complete$washu_assigned_ethnicity = NULL
clin_merge_subtype_avail_complete$Sample = NULL

clin_merge_subtype_avail_complete[ clin_merge_subtype_avail_complete == "NA" ] <- NA
clin_merge_subtype_avail_complete[ clin_merge_subtype_avail_complete == "[Not Available]"] <- NA

clin_merge_demo = clin_merge_subtype_avail_complete
clin_merge_demo = clin_merge_demo[complete.cases(clin_merge_demo),] 

clin_merge_demo$age_at_initial_pathologic_diagnosis = as.numeric(clin_merge_demo$age_at_initial_pathologic_diagnosis)

demo_table_df <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(demo_table_df) <- c("cancer","Sample_size","Female_ratio","Avg_onset_age","young_ct","reglr_ct","Young_adult_ratio")

for (cancer in unique(clin_merge_demo$acronym)) {
  temp = clin_merge_demo[clin_merge_demo$acronym==cancer,]
  
  Sample_size = nrow(temp)
  
  male_ct = sum(temp$gender=="MALE")
  fmle_ct = sum(temp$gender=="FEMALE")
  Female_ratio = round(fmle_ct/(male_ct+fmle_ct), digits=2)
  
  Avg_onset_age = round(sum(temp$age_at_initial_pathologic_diagnosis)/Sample_size, digits=1)
  
  young_ct = sum(temp$age_binary=="TRUE")
  reglr_ct = sum(temp$age_binary=="FALSE")
  
  Young_adult_ratio = round(young_ct/(reglr_ct+young_ct), digits=2)
  
  demo_table_row = cbind(cancer, Sample_size, Female_ratio, Avg_onset_age, young_ct, reglr_ct, Young_adult_ratio)
  
  demo_table_df <- rbind(demo_table_df, demo_table_row)
}

colnames(demo_table_df) <- c("Abbr.","Sample size","Female ratio","Mean onset age","Young adult cases","Later onset cases","Young adult ratio")

demo_table_df$Abbr. = as.character(demo_table_df$Abbr.)
demo_table_df = demo_table_df[order(demo_table_df$Abbr.),]

rownames(demo_table_df) <- demo_table_df$Abbr.

demo_table_df$Abbr. = NULL

demo_table_df$`Young adult cases` = as.character(demo_table_df$`Young adult cases`)
demo_table_df$`Young adult cases` = as.numeric(demo_table_df$`Young adult cases`)

demo_table_df$`Later onset cases` = as.character(demo_table_df$`Later onset cases`)
demo_table_df$`Later onset cases` = as.numeric(demo_table_df$`Later onset cases`)

demo_table_fin = demo_table_df[demo_table_df$`Young adult cases` >= 40  & demo_table_df$`Later onset cases` >= 40,]

Cancer = c("Breast invasive carcinoma","Cervical squamous cell carcinoma and endocervical adenocarcinoma",
            "Colon adenocarcinoma","Head and Neck squamous cell carcinoma","Kidney renal clear cell carcinoma",
            "Kidney renal papillary cell carcinoma","Brain Lower Grade Glioma","Liver hepatocellular carcinoma",
            "Ovarian serous cystadenocarcinoma","Pheochromocytoma and Paraganglioma","Sarcoma","Skin Cutaneous Melanoma",
            "Thyroid carcinoma","Uterine Corpus Endometrial Carcinoma")

demo_table_fin = add_column(demo_table_fin, Cancer, .before ="Sample size")

demographics = kable(demo_table_fin,"html",align="l") %>% 
  kable_styling("striped",full_width = FALSE)

demographics

############################### SUBTYPE PLOTS ###############################

library(ggplot2)
library(scales)

# BRCA (multiple)
# CESC (squamous vs adeno)
# COAD (multiple)
# HNSC (POS VS NEG)
# LGG  (multiple)
# SARC (multiple)
# UCEC (multiple)

clin_subtype_df = clin_merge_demo

clin_subtype_df$plot_age[clin_subtype_df$age_binary==TRUE] = "<= 50"
clin_subtype_df$plot_age[clin_subtype_df$age_binary==FALSE] = "> 50"

clin_subtype_female = clin_subtype_df[clin_subtype_df$acronym=="BRCA" | clin_subtype_df$acronym=="CESC" | clin_subtype_df$acronym=="UCEC",]

clin_subtype_allgen = clin_subtype_df[clin_subtype_df$acronym=="COAD" | clin_subtype_df$acronym=="HNSC" |
                                      clin_subtype_df$acronym=="LGG"  | clin_subtype_df$acronym=="SARC",]

##############################################################################################################################
for (cancer in unique(clin_subtype_female$acronym)) {
  for (subt in unique(clin_subtype_female$SUBTYPE)) {
    
    clin_subtype_female$subt_prop[clin_subtype_female$acronym==cancer & clin_subtype_female$SUBTYPE==subt & clin_subtype_female$age_binary==TRUE] =
      nrow(clin_subtype_female[clin_subtype_female$acronym==cancer & clin_subtype_female$SUBTYPE==subt & clin_subtype_female$age_binary==TRUE,])/
      nrow(clin_subtype_female[clin_subtype_female$acronym==cancer & clin_subtype_female$age_binary==TRUE,])
    
    clin_subtype_female$subt_prop[clin_subtype_female$acronym==cancer & clin_subtype_female$SUBTYPE==subt & clin_subtype_female$age_binary==FALSE] =
      nrow(clin_subtype_female[clin_subtype_female$acronym==cancer & clin_subtype_female$SUBTYPE==subt & clin_subtype_female$age_binary==FALSE,])/
      nrow(clin_subtype_female[clin_subtype_female$acronym==cancer & clin_subtype_female$age_binary==FALSE,])
  }
}

clin_subtype_female$subt_perc = round(clin_subtype_female$subt_prop*100, digits=1)

p = ggplot(data=clin_subtype_female, aes(x = plot_age, fill = SUBTYPE))
p = p + facet_wrap( .~acronym, drop=T,scale="free", nrow=3)
p = p + geom_bar(position = "fill")
p = p + scale_fill_brewer(palette = "Set3")
p = p + labs(x = "age at pathologic diagnosis (yrs.)", #y = "Subtype (%)",
             fill = "TCGA Subtype") + theme_bw()
p = p + scale_y_continuous(labels=percent)
p = p + coord_flip()
p = p + ylab(element_blank()) ###
p
fn = "out/Subtype_Female_hBarplotBinary_Percent.pdf"
ggsave(fn, w=5, h=4, useDingbats=FALSE)
##############################################################################################################################
for (cancer in unique(clin_subtype_allgen$acronym)) {
  for (subt in unique(clin_subtype_allgen$SUBTYPE)) {
    
    clin_subtype_allgen$subt_prop[clin_subtype_allgen$acronym==cancer & clin_subtype_allgen$SUBTYPE==subt & clin_subtype_allgen$age_binary==TRUE] =
      nrow(clin_subtype_allgen[clin_subtype_allgen$acronym==cancer & clin_subtype_allgen$SUBTYPE==subt & clin_subtype_allgen$age_binary==TRUE,])/
      nrow(clin_subtype_allgen[clin_subtype_allgen$acronym==cancer & clin_subtype_allgen$age_binary==TRUE,])
    
    clin_subtype_allgen$subt_prop[clin_subtype_allgen$acronym==cancer & clin_subtype_allgen$SUBTYPE==subt & clin_subtype_allgen$age_binary==FALSE] =
      nrow(clin_subtype_allgen[clin_subtype_allgen$acronym==cancer & clin_subtype_allgen$SUBTYPE==subt & clin_subtype_allgen$age_binary==FALSE,])/
      nrow(clin_subtype_allgen[clin_subtype_allgen$acronym==cancer & clin_subtype_allgen$age_binary==FALSE,])
  }
}

clin_subtype_allgen$subt_perc = round(clin_subtype_allgen$subt_prop*100, digits=1)

manualFillGeneVar = colorRampPalette(brewer.pal(12, "Set3"))(length(unique(clin_subtype_allgen$SUBTYPE))) 
p = ggplot(data=clin_subtype_allgen, aes(x = plot_age, fill = SUBTYPE))
p = p + facet_wrap( .~acronym, drop=T,scale="free", nrow=4)
p = p + geom_bar(position = "fill")
p = p + scale_fill_manual(values=manualFillGeneVar)
p = p + labs(x = "age at pathologic diagnosis (yrs.)", #y = "Subtype (%)",
             fill = "TCGA Subtype") + theme_bw()
p = p + scale_y_continuous(labels=percent)
p = p + coord_flip()
p = p + ylab(element_blank()) ###
p
fn = "out/Subtype_hBarplotBinary_Percent.pdf"
ggsave(fn, w=5, h=5, useDingbats=FALSE)

######################  FISHER'S EXACT SUBTYPES ######################

clin_merge_demo_UCEC_old_notPOLE = sum(clin_merge_demo$acronym=="UCEC" & clin_merge_demo$age_binary==FALSE & clin_merge_demo$SUBTYPE!="POLE")
clin_merge_demo_UCEC_yng_notPOLE = sum(clin_merge_demo$acronym=="UCEC" & clin_merge_demo$age_binary==TRUE & clin_merge_demo$SUBTYPE!="POLE")

clin_merge_demo_UCEC_old_POLE = sum(clin_merge_demo$acronym=="UCEC" & clin_merge_demo$age_binary==FALSE & clin_merge_demo$SUBTYPE=="POLE")
clin_merge_demo_UCEC_yng_POLE = sum(clin_merge_demo$acronym=="UCEC" & clin_merge_demo$age_binary==TRUE & clin_merge_demo$SUBTYPE=="POLE")

old_set = c(clin_merge_demo_UCEC_old_notPOLE, clin_merge_demo_UCEC_old_POLE)
yng_set = c(clin_merge_demo_UCEC_yng_notPOLE, clin_merge_demo_UCEC_yng_POLE)

UCEC_POLE_fisher = cbind(old_set, yng_set)

fisher.test(UCEC_POLE_fisher)
#------------------------------------------------------------------------------------------#
clin_merge_demo_UCEC_old_notCN_LOW = sum(clin_merge_demo$acronym=="UCEC" & clin_merge_demo$age_binary==FALSE & clin_merge_demo$SUBTYPE!="CN_LOW")
clin_merge_demo_UCEC_yng_notCN_LOW = sum(clin_merge_demo$acronym=="UCEC" & clin_merge_demo$age_binary==TRUE & clin_merge_demo$SUBTYPE!="CN_LOW")

clin_merge_demo_UCEC_old_CN_LOW = sum(clin_merge_demo$acronym=="UCEC" & clin_merge_demo$age_binary==FALSE & clin_merge_demo$SUBTYPE=="CN_LOW")
clin_merge_demo_UCEC_yng_CN_LOW = sum(clin_merge_demo$acronym=="UCEC" & clin_merge_demo$age_binary==TRUE & clin_merge_demo$SUBTYPE=="CN_LOW")

old_set = c(clin_merge_demo_UCEC_old_notCN_LOW, clin_merge_demo_UCEC_old_CN_LOW)
yng_set = c(clin_merge_demo_UCEC_yng_notCN_LOW, clin_merge_demo_UCEC_yng_CN_LOW)

UCEC_CN_LOW_fisher = cbind(old_set, yng_set)

fisher.test(UCEC_CN_LOW_fisher)
#------------------------------------------------------------------------------------------#
clin_merge_demo_UCEC_old_notCN_HIGH = sum(clin_merge_demo$acronym=="UCEC" & clin_merge_demo$age_binary==FALSE & clin_merge_demo$SUBTYPE!="CN_HIGH")
clin_merge_demo_UCEC_yng_notCN_HIGH = sum(clin_merge_demo$acronym=="UCEC" & clin_merge_demo$age_binary==TRUE & clin_merge_demo$SUBTYPE!="CN_HIGH")

clin_merge_demo_UCEC_old_CN_HIGH = sum(clin_merge_demo$acronym=="UCEC" & clin_merge_demo$age_binary==FALSE & clin_merge_demo$SUBTYPE=="CN_HIGH")
clin_merge_demo_UCEC_yng_CN_HIGH = sum(clin_merge_demo$acronym=="UCEC" & clin_merge_demo$age_binary==TRUE & clin_merge_demo$SUBTYPE=="CN_HIGH")

old_set = c(clin_merge_demo_UCEC_old_notCN_HIGH, clin_merge_demo_UCEC_old_CN_HIGH)
yng_set = c(clin_merge_demo_UCEC_yng_notCN_HIGH, clin_merge_demo_UCEC_yng_CN_HIGH)

UCEC_CN_HIGH_fisher = cbind(old_set, yng_set)

fisher.test(UCEC_CN_HIGH_fisher)
#------------------------------------------------------------------------------------------#
clin_merge_demo_BRCA_old_notBasal = sum(clin_merge_demo$acronym=="BRCA" & clin_merge_demo$age_binary==FALSE & clin_merge_demo$SUBTYPE!="Basal")
clin_merge_demo_BRCA_yng_notBasal = sum(clin_merge_demo$acronym=="BRCA" & clin_merge_demo$age_binary==TRUE & clin_merge_demo$SUBTYPE!="Basal")

clin_merge_demo_BRCA_old_Basal = sum(clin_merge_demo$acronym=="BRCA" & clin_merge_demo$age_binary==FALSE & clin_merge_demo$SUBTYPE=="Basal")
clin_merge_demo_BRCA_yng_Basal = sum(clin_merge_demo$acronym=="BRCA" & clin_merge_demo$age_binary==TRUE & clin_merge_demo$SUBTYPE=="Basal")

old_set = c(clin_merge_demo_BRCA_old_notBasal, clin_merge_demo_BRCA_old_Basal)
yng_set = c(clin_merge_demo_BRCA_yng_notBasal, clin_merge_demo_BRCA_yng_Basal)

BRCA_Basal_fisher = cbind(old_set, yng_set)

fisher.test(BRCA_Basal_fisher)
#------------------------------------------------------------------------------------------#
clin_merge_demo$SUBTYPE[clin_merge_demo$SUBTYPE=="MFS/UPS"] = "MFS_UPS"

clin_merge_demo_SARC_old_notMFS_UPS = sum(clin_merge_demo$acronym=="SARC" & clin_merge_demo$age_binary==FALSE & clin_merge_demo$SUBTYPE!="MFS_UPS")
clin_merge_demo_SARC_yng_notMFS_UPS = sum(clin_merge_demo$acronym=="SARC" & clin_merge_demo$age_binary==TRUE & clin_merge_demo$SUBTYPE!="MFS_UPS")

clin_merge_demo_SARC_old_MFS_UPS = sum(clin_merge_demo$acronym=="SARC" & clin_merge_demo$age_binary==FALSE & clin_merge_demo$SUBTYPE=="MFS_UPS")
clin_merge_demo_SARC_yng_MFS_UPS = sum(clin_merge_demo$acronym=="SARC" & clin_merge_demo$age_binary==TRUE & clin_merge_demo$SUBTYPE=="MFS_UPS")

old_set = c(clin_merge_demo_SARC_old_notMFS_UPS, clin_merge_demo_SARC_old_MFS_UPS)
yng_set = c(clin_merge_demo_SARC_yng_notMFS_UPS, clin_merge_demo_SARC_yng_MFS_UPS)

SARC_MFS_UPS_fisher = cbind(old_set, yng_set)

fisher.test(SARC_MFS_UPS_fisher)
#------------------------------------------------------------------------------------------#
clin_merge_demo$SUBTYPE[clin_merge_demo$SUBTYPE=="IDHmut-non-codel"] = "IDHmut_non_codel"

clin_merge_demo_LGG_old_notIDHmut_non_codel = sum(clin_merge_demo$acronym=="LGG" & clin_merge_demo$age_binary==FALSE & clin_merge_demo$SUBTYPE!="IDHmut_non_codel")
clin_merge_demo_LGG_yng_notIDHmut_non_codel = sum(clin_merge_demo$acronym=="LGG" & clin_merge_demo$age_binary==TRUE & clin_merge_demo$SUBTYPE!="IDHmut_non_codel")

clin_merge_demo_LGG_old_IDHmut_non_codel = sum(clin_merge_demo$acronym=="LGG" & clin_merge_demo$age_binary==FALSE & clin_merge_demo$SUBTYPE=="IDHmut_non_codel")
clin_merge_demo_LGG_yng_IDHmut_non_codel = sum(clin_merge_demo$acronym=="LGG" & clin_merge_demo$age_binary==TRUE & clin_merge_demo$SUBTYPE=="IDHmut_non_codel")

old_set = c(clin_merge_demo_LGG_old_notIDHmut_non_codel, clin_merge_demo_LGG_old_IDHmut_non_codel)
yng_set = c(clin_merge_demo_LGG_yng_notIDHmut_non_codel, clin_merge_demo_LGG_yng_IDHmut_non_codel)

LGG_IDHmut_non_codel_fisher = cbind(old_set, yng_set)

fisher.test(LGG_IDHmut_non_codel_fisher)
#------------------------------------------------------------------------------------------#
clin_merge_demo_LGG_old_notIDHwt = sum(clin_merge_demo$acronym=="LGG" & clin_merge_demo$age_binary==FALSE & clin_merge_demo$SUBTYPE!="IDHwt")
clin_merge_demo_LGG_yng_notIDHwt = sum(clin_merge_demo$acronym=="LGG" & clin_merge_demo$age_binary==TRUE & clin_merge_demo$SUBTYPE!="IDHwt")

clin_merge_demo_LGG_old_IDHwt = sum(clin_merge_demo$acronym=="LGG" & clin_merge_demo$age_binary==FALSE & clin_merge_demo$SUBTYPE=="IDHwt")
clin_merge_demo_LGG_yng_IDHwt = sum(clin_merge_demo$acronym=="LGG" & clin_merge_demo$age_binary==TRUE & clin_merge_demo$SUBTYPE=="IDHwt")

old_set = c(clin_merge_demo_LGG_old_notIDHwt, clin_merge_demo_LGG_old_IDHwt)
yng_set = c(clin_merge_demo_LGG_yng_notIDHwt, clin_merge_demo_LGG_yng_IDHwt)

LGG_IDHwt_fisher = cbind(old_set, yng_set)

fisher.test(LGG_IDHwt_fisher)
#------------------------------------------------------------------------------------------#
# clin_merge_demo$SUBTYPE[clin_merge_demo$SUBTYPE=="IDHmut-codel"] = "IDHmut_codel" ## not significant!
# data:  LGG_IDHmut_codel_fisher
# p-value = 0.06444

######################  CNV  ###################### 

# source the files with plotting related functions
source("../global_aes_out.R")

cnv_seg = "TCGA_mastercalls.abs_segtabs.fixed.CNVperSample.tsv"
CNV_data = read.table(header=T, quote = "", sep="\t", file = cnv_seg, stringsAsFactors=FALSE)
colnames(CNV_data) = c("CNV_segs","SAMPLE_BARCODE")

CNV_seg_merge_df = merge(clin_merge_demo,CNV_data, by="SAMPLE_BARCODE")

CNV_seg_merge_df$plot_age[CNV_seg_merge_df$age_binary==TRUE] = "<= 50"
CNV_seg_merge_df$plot_age[CNV_seg_merge_df$age_binary==FALSE] = "> 50"

CNV_seg_merge_df = CNV_seg_merge_df[CNV_seg_merge_df$acronym == "BRCA" | CNV_seg_merge_df$acronym == "CESC" |
                                    CNV_seg_merge_df$acronym == "COAD" | CNV_seg_merge_df$acronym == "HNSC" |
                                    CNV_seg_merge_df$acronym == "KIRC" | CNV_seg_merge_df$acronym == "KIRP" |
                                    CNV_seg_merge_df$acronym == "LGG"  | CNV_seg_merge_df$acronym == "LIHC" |
                                    CNV_seg_merge_df$acronym == "OV"   | CNV_seg_merge_df$acronym == "PCPG" |
                                    CNV_seg_merge_df$acronym == "SARC" | CNV_seg_merge_df$acronym == "SKCM" | 
                                    CNV_seg_merge_df$acronym == "THCA" | CNV_seg_merge_df$acronym == "UCEC",]

p1 = ggplot(data=CNV_seg_merge_df,aes(x=plot_age,y=log10(CNV_segs),fill=plot_age))
p1 = p1 + geom_violin(trim = TRUE)
p1 = p1 + facet_grid( .~acronym, drop=T, space ="free", scale = "free")
#p1 = p1 + geom_jitter(aes(color = plot_age), size=0.35, alpha=1)
p1 = p1 + scale_x_discrete(drop=FALSE) #+ scale_y_discrete(drop=FALSE)
p1 = p1 + ylab("log10(CNV segments)") + xlab("age at pathologic diagnosis (yrs.)")
#p1 = p1 + stat_summary(fun.data=mean_sdl ,geom="pointrange", color="black")
p1 = p1 + theme_bw()
p1 = p1 + scale_fill_manual(values=c("<= 50" = "#00BFC4", "> 50" = "#F8766D"))
p1 = p1 + theme(legend.position="none")
p1 = p1 + stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.75)
p1
fn = 'out/YoungAdult_CNV_Segment_Violin.pdf'
ggsave(fn,h=3.5, w=11.5,useDingbats=F)

########################################################################################

p = ggplot(data=CNV_seg_merge_df)
p = p + stat_summary(aes(x = plot_age, y = CNV_segs), geom = "col", fun = sum,
                     colour = "black", fill = "dodgerblue2")
p = p + facet_wrap( .~acronym, drop=T,scale="fixed", nrow=1)
#p = p + geom_col()#stat = "count") #position = "stack")
p = p + labs(y = "CNV segments", x = "age at pathologic diagnosis (yrs.)") + theme_bw()
p = p + theme(legend.position="top")
#p = p + scale_color_manual(values=c("<= 50" = "#00BFC4", "> 50" = "#F8766D"))
#p = p + theme(legend.title = element_text(size=8))#, legend.text = element_text(size=8))
#p = p + xlab(element_blank()) ###
p
fn = "out/YoungAdult_CNV_Segment_Barplot.pdf"
ggsave(fn, w=11, h=2.5, useDingbats=FALSE)




# Code that might be useful later #
# clin_merge_demo = clin_merge_subtype_avail_complete[,c("bcr_patient_barcode","acronym","age_at_initial_pathologic_diagnosis",
#                                                        "gender","SUBTYPE","age_binary")]

#demo_table_row = as.data.frame(cbind(cancer, Sample_size, Female_ratio, Avg_onset_age, young_ct, reglr_ct, Young_adult_ratio))
#demographics = formattable(demo_table_df)

# BRCA_subtype = clin_subtype_df[clin_subtype_df$acronym=="BRCA",]
# CESC_subtype = clin_subtype_df[clin_subtype_df$acronym=="CESC",]

# # SUBTYPE COUNT PLOT #
# p1 = ggplot(data=clin_subtype_df,aes(x=plot_age,y=SUBTYPE))#,colour = age_binary))
# p1 = p1 + geom_count()
# p1 = p1 + facet_wrap(.~acronym, drop=T, scale = "free", nrow = 7)
# p1 = p1 + scale_x_discrete(drop=FALSE) #+ scale_y_discrete(drop=FALSE)
# p1 = p1 + ylab("TCGA subtype") + xlab("age at pathologic diagnosis (yrs.)")
# p1 = p1 + theme_bw()
# #p1 = p1 + scale_color_manual(values=c("young_skewed" = "#00BFC4", "regular_skewed" = "#F8766D"))
# p1 = p1 + theme(legend.position="none")
# p1
# fn = 'out/Subtype_OnsetBinaryCount.pdf'
# ggsave(fn,h=12, w=5,useDingbats=F)

#p = p + geom_text(aes(label=subt_perc, y=subt_prop))

#   