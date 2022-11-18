### phosphorylation abundance distribution
rm(list=ls())
library(ggplot2)
setwd('/home/shengli/projects/AEG_proteomics/data/phosphoproteome')
prot_ibaq <- read.table('Phosphoproteomics_iBAQ103_log2quantile_normlization_df.txt',header=T)
sample_cancer <- as.character(unique(prot_ibaq[which(prot_ibaq[,'Group']=='Cancer'),'Sample']))
sample_normal <- as.character(unique(prot_ibaq[which(prot_ibaq[,'Group']=='Normal'),'Sample']))
samples <- c(sample_cancer,sample_normal)
pdf('/home/shengli/projects/AEG_proteomics/figures/Phosphoproteomics_iBAQ_norm_distribution.pdf',width=20,height=5)
ggplot(prot_ibaq,aes(x=Sample,y=Abundance,fill=Group)) +
  geom_boxplot()
dev.off()
