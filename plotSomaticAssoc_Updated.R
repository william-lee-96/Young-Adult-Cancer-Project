##### plotSomaticAssoc_Updated.R #####
# Updated by Will Lee, August 2020

bdir = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/somatic_mutation"
setwd(bdir)

# source the files with plotting related functions
source("../global_aes_out.R")

### associations - updated by Will Lee ###
assoc_somatic_n = "out/somatic_likelyfunctional_driver_onsetAgeBinary_assoc_Updated.txt"
assoc_somatic = read.table(sep="\t",header=T,file=assoc_somatic_n , stringsAsFactors=FALSE)

### Will's updated glm plots ###

# bubble plot, somatic #
p = ggplot(data=assoc_somatic,aes(x=cancer,y = coefficient,color = cancer))
p = p + geom_point(aes(size = -log10(FDR)),alpha=0.8)
p = p + geom_text_repel(aes(label=ifelse(FDR<0.05,gene,NA)), size=3, colour="black", show.legend = FALSE)
p = p + getPCACancerColor()
p = p + labs(x="Cancer",y= "Coefficient") + ylim(-2.5,2.6)
p = p + theme(legend.position = "top")
p = p + geom_hline(yintercept = 0, alpha=0.1)
p = p + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12, angle= 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())
p
fn = 'out/SomaticAssocOnsetBinaryBubble_Updated.pdf'
ggsave(fn,h=5.75, w=5,useDingbat=F)