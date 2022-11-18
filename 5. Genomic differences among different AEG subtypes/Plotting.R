### plot mutation effect on protein expression
library(ggplot2)
setwd('/home/shengli/projects/AEG_proteomics/results/WES/2022.08')
mut_prot <- read.table('AEG_specific_mut_affect_proteins_sig.txt',header=T,sep='\t')
mut_prot$lgfdr <- -log10(mut_prot[,'FDR'])
mut_prot[which(mut_prot[,'lgfdr'] > 15),'lgfdr'] <- 15

s1_mutgenes <- c('LEPR','SORCS3','POTEE','COL2A1','DLG2')
s1_proteins <- c('VHL','YIF1B','MAZ','PHF21A','CRCP',
                 'POPDC2','TMEM35B','PTP4A3','GLT8D2','CYS1',
                 'FRA10AC1','LRIG3','FABP6','KLK6','TPX2',
                 'LRRC32','ARL5B','CEACAM7','ACTR6',
                 'RHOBTB3','NAT8','GPS2','BAZ2A','TMEM135')
s2_mutgenes <- c('NCKAP1','DOPEY1','DENND1B','TRPC4','MAP7')
s2_proteins <- c('CHRDL2','PTGS2','IFT22','TASP1','C1QTNF5',
                 'L1RE1','LIMS4','PLAC8','CCR1',
                 'H1-1','SFRP2','SH2D2A',
                 'CLRN3','RBM38','HLA-DOA','CALHM6','HERPUD2',
                 'GLYCTK','PHLDA2','SFRP4','PFKFB3')
s3_mutgenes <- c('WIZ','LGALS9','CALCRL','FGF18','REEP6')
s3_proteins <- c('CXCL8','VHL','CPXM1','RBP2','SFRP4',
                 'SLC15A4','CEBPG','PLA1A','CKS2','SLC2A5',
                 'CCDC115','RAB11FIP3','LHFPL2','OLR1',
                 'RAMAC','UGT2B7','DRP2','TM4SF20',
                 'CLDN6','CWF19L2','B4GALT7','CCNC','IDS')

sel_s3 <- mut_prot[which((mut_prot[,'Gene'] %in% s3_mutgenes & mut_prot[,'Protein_name'] %in% s3_proteins) & mut_prot[,'FDR'] < 0.05),]

pdf('/home/shengli/projects/AEG_proteomics/figures/2022.08/S-III_mut_protein_top5.pdf',width=10,height=3.5)
ggplot(sel_s3,aes(x=Protein_name,y=Gene))+
  geom_point(aes(size=lgfdr,col=FoldChange))+
  scale_color_gradient2(low="blue",high="red",na.value="white",name="Fold change difference")+
  scale_y_discrete(limit=s3_mutgenes,expand=c(0.05,0.05))+
  scale_x_discrete(limit=s3_proteins,expand=c(0.05,0.05)) +
  theme(panel.background=element_rect(colour="black",fill="white"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=11,colour="black"),
        axis.text.x=element_text(size=10,colour="black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.key=element_rect(fill="white",colour="black"))
dev.off()

