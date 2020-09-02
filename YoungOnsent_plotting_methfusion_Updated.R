##### YoungOnsent_plotting_methfusion_Updated.R #####
# Updated by Will Lee, August 2020

bdir = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/methylation_fusion"
setwd(bdir)

# source the files with plotting related functions
source("../global_aes_out.R")

### associations - updated by Will Lee ###
assoc_meth_n = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/methylation_fusion/out/meth_events_onsetAgeBinary_assoc_Updated.txt"
assoc_meth = read.table(sep="\t",header=T,file=assoc_meth_n , stringsAsFactors=FALSE)

assoc_fusion_n = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/methylation_fusion/out/fus_events_onsetAgeBinary_assoc_Updated.txt"
assoc_fusion = read.table(sep="\t",header=T,file=assoc_fusion_n , stringsAsFactors=FALSE)

assoc_cnv_n = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/methylation_fusion/out/cnv_events_onsetAgeBinary_assoc_Updated.txt"
assoc_cnv = read.table(sep="\t",header=T,file=assoc_cnv_n , stringsAsFactors=FALSE)

############################################################################################

cnv_log_adj = assoc_cnv

cnv_log_adj$neg_vals = cnv_log_adj$coefficient < 0
cnv_log_adj$coefficient = abs(cnv_log_adj$coefficient)
cnv_log_adj$coefficient = log10(cnv_log_adj$coefficient)

cnv_log_adj$coefficient[cnv_log_adj$coefficient > 0 & cnv_log_adj$neg_vals==TRUE] = 
cnv_log_adj$coefficient[cnv_log_adj$coefficient > 0 & cnv_log_adj$neg_vals==TRUE]*(-1)

cnv_log_adj$coefficient[cnv_log_adj$coefficient < 0 & cnv_log_adj$neg_vals==FALSE] = 
cnv_log_adj$coefficient[cnv_log_adj$coefficient < 0 & cnv_log_adj$neg_vals==FALSE]*(-1)

############################################################################################

meth_log_adj = assoc_meth

meth_log_adj$neg_vals = meth_log_adj$coefficient < 0
meth_log_adj$coefficient = abs(meth_log_adj$coefficient)
meth_log_adj$coefficient = log10(meth_log_adj$coefficient)

meth_log_adj$coefficient[meth_log_adj$coefficient > 0 & meth_log_adj$neg_vals==TRUE] = 
meth_log_adj$coefficient[meth_log_adj$coefficient > 0 & meth_log_adj$neg_vals==TRUE]*(-1)

meth_log_adj$coefficient[meth_log_adj$coefficient < 0 & meth_log_adj$neg_vals==FALSE] = 
meth_log_adj$coefficient[meth_log_adj$coefficient < 0 & meth_log_adj$neg_vals==FALSE]*(-1)

############################################################################################

fusion_log_adj = assoc_fusion

fusion_log_adj$neg_vals = fusion_log_adj$coefficient < 0
fusion_log_adj$coefficient = abs(fusion_log_adj$coefficient)
fusion_log_adj$coefficient = log10(fusion_log_adj$coefficient)

fusion_log_adj$coefficient[fusion_log_adj$coefficient > 0 & fusion_log_adj$neg_vals==TRUE] = 
fusion_log_adj$coefficient[fusion_log_adj$coefficient > 0 & fusion_log_adj$neg_vals==TRUE]*(-1)

fusion_log_adj$coefficient[fusion_log_adj$coefficient < 0 & fusion_log_adj$neg_vals==FALSE] = 
fusion_log_adj$coefficient[fusion_log_adj$coefficient < 0 & fusion_log_adj$neg_vals==FALSE]*(-1)

############################################################################################

### Will's updated glm plots ###

# bubble plot, methylation #
p = ggplot(data=meth_log_adj,aes(x=cancer,y = coefficient,color = cancer))
p = p + geom_point(aes(size = -log10(FDR)),alpha=0.8)
p = p + geom_text_repel(aes(label=ifelse(FDR<0.15,meth_event,NA)), size=2.75, color="black",show.legend = FALSE)
p = p + getPCACancerColor()
p = p + labs(x="Cancer",y= "log10(Coefficient)") + ylim(-1.5,0.75)
p = p + geom_hline(yintercept = 0, alpha=0.1)
p = p  + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12, angle= 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
fn = 'out/MethAssocOnsetBinaryBubble_Updated.pdf'
ggsave(fn,h=5.75, w=5,useDingbat=F)

# bubble plot, fusion #
p = ggplot(data=fusion_log_adj,aes(x=cancer,y = coefficient,color = cancer))
p = p + geom_point(aes(size = -log10(FDR)),alpha=0.8)
p = p + geom_text_repel(aes(label=ifelse(FDR<0.15,fus_event,NA)), size=3, colour="black", show.legend = FALSE)
p = p + getPCACancerColor()
p = p + labs(x="Cancer",y= "log10(Coefficient)") + ylim(-0.5,1.5)
p = p + geom_hline(yintercept = 0, alpha=0.1)
p = p  + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12, angle= 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
fn = 'out/FusionAssocOnsetBinaryBubble_Updated.pdf'
ggsave(fn,h=5.75, w=3.5,useDingbat=F)

# bubble plot, cnv #
p = ggplot(data=cnv_log_adj,aes(x=cancer,y = coefficient,color = cancer))
p = p + geom_point(aes(size = -log10(FDR)),alpha=0.8)
p = p + geom_text_repel(aes(label=ifelse(FDR<0.15,cnv_event,NA)), size=2.75, colour="black", show.legend = FALSE)
p = p + getPCACancerColor()
p = p + labs(x="Cancer",y= "log10(Coefficient)") + ylim(-1.5,1.5)
p = p + geom_hline(yintercept = 0, alpha=0.1)
p = p  + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12, angle= 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
fn = 'out/CNVAssocOnsetBinaryBubble_Updated.pdf'
ggsave(fn,h=5.75, w=5,useDingbat=F)