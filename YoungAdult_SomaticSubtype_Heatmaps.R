# William Lee @ May 2020
# work in progress

bdir = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/somatic_mutation"
setwd(bdir)

# source the files with plotting related functions
source("../global_aes_out.R")

### MAIN ###

### MC3 mutation call file (only include the likely driver for the first-ass analysis) ###
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

##### statistical testing: pan-gene analysis within each cancer type #####
# somatic mutation in driver genes ~ onset age + covariates ("gender","Subtype","PC1","PC2")

clin_merge_subtype_avail_complete$cancer = NULL
clin_merge_subtype_avail_complete$ethnicity = NULL
clin_merge_subtype_avail_complete$washu_assigned_ethnicity = NULL
clin_merge_subtype_avail_complete$Sample = NULL
clin_merge_subtype_avail_complete[ clin_merge_subtype_avail_complete == "NA" ] <- NA
clin_merge_subtype_avail_complete[ clin_merge_subtype_avail_complete == "[Not Available]"] <- NA

clin_somatic_mut_ggplot <- data.frame(matrix(ncol = 35, nrow = 0))
colnames(clin_somatic_mut_ggplot) <- c("bcr_patient_barcode", "acronym", "age_at_initial_pathologic_diagnosis", "gender",
                                       "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12",
                                       "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "SAMPLE_BARCODE", "DISEASE",
                                       "SUBTYPE", "age_binary", "Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", "HGVSp_Short")

colnames(clin_somatic_mut_ggplot)[33] = "plot_age"
colnames(clin_somatic_mut_ggplot)[34] = "gene"
colnames(clin_somatic_mut_ggplot)[35] = "gene_m_stat"

# cancer = "BRCA" ; gene = "TP53" # can use to trouble-shoot
# conduct tests by each cancer (denoted as acronym)
for (cancer in unique(clin_merge_subtype_avail_complete$acronym)){
  
  #if (cancer == "BRCA" | cancer == "CESC" | cancer == "COAD" | cancer == "LGG" | cancer == "HNSC"| cancer == "SKCM" | cancer == "UCEC") {
  if (cancer == "BRCA") {
    
    clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete[clin_merge_subtype_avail_complete$acronym==cancer,]
    clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete_c[complete.cases(clin_merge_subtype_avail_complete_c),]
    
    # conduct the test by gene
    for (gene in unique(somatic_likelyfunctional_driver$Hugo_Symbol)){
      
      # subset mutation to that gene and get unique mutation
      somatic_likelyfunctional_driver_g = somatic_likelyfunctional_driver[somatic_likelyfunctional_driver$Hugo_Symbol==gene,]
      somatic_likelyfunctional_driver_g_uniq = somatic_likelyfunctional_driver_g[!duplicated(somatic_likelyfunctional_driver_g$bcr_patient_barcode),]
      clin_merge_subtype_avail_complete_c_merge_somatic = merge(clin_merge_subtype_avail_complete_c,somatic_likelyfunctional_driver_g_uniq, by= "bcr_patient_barcode", all.x=T, all.y = F)
      
      # model whether the sample carry somatic mutation in the gene (Hugo Symbol is a standardized gene name)
      # samples with mutations considered as 1, no mutations considered as 0
      clin_merge_subtype_avail_complete_c_merge_somatic$Hugo_Symbol[!is.na(clin_merge_subtype_avail_complete_c_merge_somatic$Hugo_Symbol)] = "MUT"
      clin_merge_subtype_avail_complete_c_merge_somatic$Hugo_Symbol[is.na(clin_merge_subtype_avail_complete_c_merge_somatic$Hugo_Symbol)] = "WT"
      
      clin_merge_subtype_avail_complete_c_merge_somatic$age_at_initial_pathologic_diagnosis = as.numeric(clin_merge_subtype_avail_complete_c_merge_somatic$age_at_initial_pathologic_diagnosis)
      
      clin_merge_subtype_avail_complete_c_merge_somatic$plot_age[clin_merge_subtype_avail_complete_c_merge_somatic$age_binary==TRUE] = "<= 50"
      clin_merge_subtype_avail_complete_c_merge_somatic$plot_age[clin_merge_subtype_avail_complete_c_merge_somatic$age_binary==FALSE] = "> 50"
      
      clin_merge_subtype_avail_complete_c_merge_somatic$gene = gene
      
      clin_merge_subtype_avail_complete_c_merge_somatic$gene_m_stat = paste(clin_merge_subtype_avail_complete_c_merge_somatic$gene,
                                                                            clin_merge_subtype_avail_complete_c_merge_somatic$Hugo_Symbol,sep=":\n")
      
      clin_somatic_mut_ggplot <- rbind(clin_somatic_mut_ggplot, clin_merge_subtype_avail_complete_c_merge_somatic)
    }
  }
}

clin_somatic_mut_ggplot_UCEC = clin_somatic_mut_ggplot[clin_somatic_mut_ggplot$acronym=="UCEC",]
clin_somatic_mut_ggplot_BRCA = clin_somatic_mut_ggplot[clin_somatic_mut_ggplot$acronym=="BRCA",]
clin_somatic_mut_ggplot_LGG = clin_somatic_mut_ggplot[clin_somatic_mut_ggplot$acronym=="LGG",]

################################# UCEC HEATMAP #################################

UCEC_hmap = clin_somatic_mut_ggplot_UCEC
total_yng = nrow(UCEC_hmap[!duplicated(UCEC_hmap$bcr_patient_barcode) & UCEC_hmap$age_binary==TRUE,])
total_late = nrow(UCEC_hmap[!duplicated(UCEC_hmap$bcr_patient_barcode) & UCEC_hmap$age_binary==FALSE,])

UCEC_hmap$hmap_label[UCEC_hmap$age_binary==TRUE] = paste("young adult, n = ", total_yng, sep="")
UCEC_hmap$hmap_label[UCEC_hmap$age_binary==FALSE] = paste("later onset, n = ", total_late, sep="")

sig_gene_mut = c("ATRX:\nMUT","BRD7:\nMUT","CNBD1:\nMUT","CTNNB1:\nMUT","FLT3:\nMUT","LATS1:\nMUT","PTEN:\nMUT","SIN3A:\nMUT")

for (gmut in sig_gene_mut) {
  for (subt in unique(UCEC_hmap$SUBTYPE)) {
    
    UCEC_hmap$hmap_data[UCEC_hmap$gene_m_stat==gmut & UCEC_hmap$SUBTYPE==subt & UCEC_hmap$age_binary==TRUE] =
      nrow(UCEC_hmap[UCEC_hmap$gene_m_stat==gmut & UCEC_hmap$SUBTYPE==subt & UCEC_hmap$age_binary==TRUE,])/total_yng
    
    UCEC_hmap$hmap_data[UCEC_hmap$gene_m_stat==gmut & UCEC_hmap$SUBTYPE==subt & UCEC_hmap$age_binary==FALSE] =
      nrow(UCEC_hmap[UCEC_hmap$gene_m_stat==gmut & UCEC_hmap$SUBTYPE==subt & UCEC_hmap$age_binary==FALSE,])/total_late
    
  }
}

UCEC_hmap = UCEC_hmap[UCEC_hmap$gene_m_stat=="ATRX:\nMUT" | UCEC_hmap$gene_m_stat=="BRD7:\nMUT" |
                        UCEC_hmap$gene_m_stat=="CNBD1:\nMUT" | UCEC_hmap$gene_m_stat=="CTNNB1:\nMUT" |
                        UCEC_hmap$gene_m_stat=="FLT3:\nMUT" | UCEC_hmap$gene_m_stat=="LATS1:\nMUT" |
                        UCEC_hmap$gene_m_stat=="PTEN:\nMUT" | UCEC_hmap$gene_m_stat=="SIN3A:\nMUT",]

UCEC_hmap$subt_perc = round(UCEC_hmap$hmap_data*100, digits=1)

# UCEC heatmap #
p = ggplot(data=UCEC_hmap, mapping = aes(x = SUBTYPE, y = gene_m_stat, fill = subt_perc))
p = p + facet_wrap( .~hmap_label, drop=T,scale="fixed", nrow=1)
p = p + geom_tile()
p = p + geom_text(aes(label=subt_perc),color="orangered",size=3.5)
p = p + theme_bw()
p = p + scale_fill_distiller(direction = 1)
p = p + theme(panel.grid.major = element_blank())
p = p + xlab(element_blank()) + ylab(element_blank())
p = p + theme(axis.text.x = element_text(angle=22.5,vjust=0.65,size=8), axis.text.y = element_text(size=8))
p = p + guides(fill=guide_colorbar(title="% UCEC"))#,barheight=7))
#p = p + theme(legend.position = "top")
p
fn = 'out/UCEC_mut_subtype_heatmap.pdf'
ggsave(fn,h=3, w=5.25,useDingbats=F)

################################# BRCA HEATMAP #################################

BRCA_hmap = clin_somatic_mut_ggplot_BRCA
total_yng = nrow(BRCA_hmap[!duplicated(BRCA_hmap$bcr_patient_barcode) & BRCA_hmap$age_binary==TRUE,])
total_late = nrow(BRCA_hmap[!duplicated(BRCA_hmap$bcr_patient_barcode) & BRCA_hmap$age_binary==FALSE,])

BRCA_hmap$hmap_label[BRCA_hmap$age_binary==TRUE] = paste("young adult, n = ", total_yng, sep="")
BRCA_hmap$hmap_label[BRCA_hmap$age_binary==FALSE] = paste("later onset, n = ", total_late, sep="")

sig_gene_mut = c("CDH1:\nMUT","GATA3:\nMUT","KMT2C:\nMUT")

for (gmut in sig_gene_mut) {
  for (subt in unique(BRCA_hmap$SUBTYPE)) {
    
    BRCA_hmap$hmap_data[BRCA_hmap$gene_m_stat==gmut & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==TRUE] =
      nrow(BRCA_hmap[BRCA_hmap$gene_m_stat==gmut & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==TRUE,])/total_yng
    
    BRCA_hmap$hmap_data[BRCA_hmap$gene_m_stat==gmut & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==FALSE] =
      nrow(BRCA_hmap[BRCA_hmap$gene_m_stat==gmut & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==FALSE,])/total_late
    
  }
}

BRCA_hmap = BRCA_hmap[BRCA_hmap$gene_m_stat=="CDH1:\nMUT" | BRCA_hmap$gene_m_stat=="GATA3:\nMUT" |
                        BRCA_hmap$gene_m_stat=="KMT2C:\nMUT",]

BRCA_hmap$subt_perc = round(BRCA_hmap$hmap_data*100, digits=1)

# BRCA heatmap #
p = ggplot(data=BRCA_hmap, mapping = aes(x = SUBTYPE, y = gene_m_stat, fill = subt_perc))
p = p + facet_wrap( .~hmap_label, drop=T,scale="fixed", nrow=1)
p = p + geom_tile()
p = p + geom_text(aes(label=subt_perc),color="orangered",size=3.5)
p = p + theme_bw()
p = p + scale_fill_distiller(direction = 1)
p = p + theme(panel.grid.major = element_blank())
p = p + xlab(element_blank()) + ylab(element_blank())
p = p + theme(axis.text.x = element_text(angle=22.5,vjust=0.65))#, axis.text.y = element_text(size=8))
p = p + guides(fill=guide_colorbar(title="% BRCA"))
#p = p + theme(legend.position = "top")
p
fn = 'out/BRCA_mut_subtype_heatmap.pdf'
ggsave(fn,h=1.75, w=5.25,useDingbats=F)

################################# LGG HEATMAP #################################

LGG_hmap = clin_somatic_mut_ggplot_LGG
total_yng = nrow(LGG_hmap[!duplicated(LGG_hmap$bcr_patient_barcode) & LGG_hmap$age_binary==TRUE,])
total_late = nrow(LGG_hmap[!duplicated(LGG_hmap$bcr_patient_barcode) & LGG_hmap$age_binary==FALSE,])

LGG_hmap$hmap_label[LGG_hmap$age_binary==TRUE] = paste("young adult, n = ", total_yng, sep="")
LGG_hmap$hmap_label[LGG_hmap$age_binary==FALSE] = paste("later onset, n = ", total_late, sep="")

sig_gene_mut = c("ATRX:\nMUT","EGFR:\nMUT","IDH1:\nMUT","TP53:\nMUT")

for (gmut in sig_gene_mut) {
  for (subt in unique(LGG_hmap$SUBTYPE)) {
    
    LGG_hmap$hmap_data[LGG_hmap$gene_m_stat==gmut & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==TRUE] =
      nrow(LGG_hmap[LGG_hmap$gene_m_stat==gmut & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==TRUE,])/total_yng
    
    LGG_hmap$hmap_data[LGG_hmap$gene_m_stat==gmut & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==FALSE] =
      nrow(LGG_hmap[LGG_hmap$gene_m_stat==gmut & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==FALSE,])/total_late
    
  }
}

LGG_hmap = LGG_hmap[LGG_hmap$gene_m_stat=="ATRX:\nMUT" | LGG_hmap$gene_m_stat=="EGFR:\nMUT" |
                      LGG_hmap$gene_m_stat=="IDH1:\nMUT" | LGG_hmap$gene_m_stat=="TP53:\nMUT",]

LGG_hmap$SUBTYPE[LGG_hmap$SUBTYPE=="IDHmut-non-codel"] = "IDHmut-non\n-codel"
LGG_hmap$SUBTYPE[LGG_hmap$SUBTYPE=="IDHmut-codel"] = "IDHmut\n-codel"

LGG_hmap$subt_perc = round(LGG_hmap$hmap_data*100, digits=1)

# LGG heatmap #
p = ggplot(data=LGG_hmap, mapping = aes(x = SUBTYPE, y = gene_m_stat, fill = subt_perc))
p = p + facet_wrap( .~hmap_label, drop=T,scale="fixed", nrow=1)
p = p + geom_tile()
p = p + geom_text(aes(label=subt_perc),color="orangered",size=3.5)
p = p + theme_bw()
p = p + scale_fill_distiller(direction = 1)
p = p + theme(panel.grid.major = element_blank())
p = p + xlab(element_blank()) + ylab(element_blank())
p = p + theme(axis.text.x = element_text(size=8))#, axis.text.y = element_text(size=8))
p = p + guides(fill=guide_colorbar(title="% LGG"))
#p = p + theme(legend.position = "top")
p
fn = 'out/LGG_mut_subtype_heatmap.pdf'
ggsave(fn,h=2.5, w=5.25,useDingbats=F)


################################# BRCA COUNT HEATMAP #################################

BRCA_hmap = clin_somatic_mut_ggplot_BRCA
total_yng = nrow(BRCA_hmap[!duplicated(BRCA_hmap$bcr_patient_barcode) & BRCA_hmap$age_binary==TRUE,])
total_late = nrow(BRCA_hmap[!duplicated(BRCA_hmap$bcr_patient_barcode) & BRCA_hmap$age_binary==FALSE,])

BRCA_hmap$hmap_label[BRCA_hmap$age_binary==TRUE] = paste("young adult, n = ", total_yng, sep="")
BRCA_hmap$hmap_label[BRCA_hmap$age_binary==FALSE] = paste("later onset, n = ", total_late, sep="")

sig_gene_mut = c("CDH1:\nMUT","GATA3:\nMUT","KMT2C:\nMUT")

for (gmut in sig_gene_mut) {
  for (subt in unique(BRCA_hmap$SUBTYPE)) {
    
    BRCA_hmap$hmap_data[BRCA_hmap$gene_m_stat==gmut & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==TRUE] =
      nrow(BRCA_hmap[BRCA_hmap$gene_m_stat==gmut & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==TRUE,])
    
    BRCA_hmap$hmap_data[BRCA_hmap$gene_m_stat==gmut & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==FALSE] =
      nrow(BRCA_hmap[BRCA_hmap$gene_m_stat==gmut & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==FALSE,])
    
  }
}

BRCA_hmap = BRCA_hmap[BRCA_hmap$gene_m_stat=="CDH1:\nMUT" | BRCA_hmap$gene_m_stat=="GATA3:\nMUT" |
                        BRCA_hmap$gene_m_stat=="KMT2C:\nMUT",]

# BRCA heatmap #
p = ggplot(data=BRCA_hmap, mapping = aes(x = SUBTYPE, y = gene_m_stat, fill = hmap_data))
p = p + facet_wrap( .~hmap_label, drop=T,scale="fixed", nrow=1)
p = p + geom_tile()
p = p + geom_text(aes(label=hmap_data),color="orangered",size=3.5)
p = p + theme_bw()
p = p + scale_fill_distiller(direction = 1)
p = p + theme(panel.grid.major = element_blank())
p = p + xlab(element_blank()) + ylab(element_blank())
p = p + theme(axis.text.x = element_text(angle=22.5,vjust=0.65))#, axis.text.y = element_text(size=8))
p = p + guides(fill=guide_colorbar(title="# BRCA"))
#p = p + theme(legend.position = "top")
p
fn = 'out/testing.pdf'
ggsave(fn,h=1.75, w=5.25,useDingbats=F)

