### waterfall plot of top mutation cancer genes, Fig. 1b
library(maftools) # version 2.6.05
setwd('/home/shengli/projects/AEG_proteomics/data/WES')
aeg <- read.maf(maf='AEG_mutation_cancer_gene.maf')
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/mutation_waterfall_cancer_genes_top30.pdf',height=6,width=10)
oncoplot(maf=aeg,draw_titv=T,top=30,showTumorSampleBarcodes=F)
dev.off()

### The number distribution of identified proteins in each sample, Fig. 1c
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
prt_num <- read.table('Proteomics_sample_numbers.txt',header=T,sep='\t')
num_tumor <- prt_num[which(prt_num[,'Type']=='tumor'),]
rownames(num_tumor) <- num_tumor[,'Sample']
num_tumor_ordered <- num_tumor[order(num_tumor[,'Protein_number'],decreasing = F),]
pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics/Protein_numbers_tumor.pdf',height=5,width=5)
plot(num_tumor_ordered[,'Protein_number'],ylim=c(8000,9500),las=1)
dev.off()
num_normal <- prt_num[which(prt_num[,'Type']=='normal'),]
rownames(num_normal) <- num_normal[,'Sample']
num_normal_ordered <- num_normal[rownames(num_tumor_ordered),]
wilcox.test(as.numeric(num_normal_ordered[,'Protein_number']),as.numeric(num_tumor_ordered[,'Protein_number']))

pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics/Protein_numbers_normal.pdf',height=5,width=5)
plot(num_normal_ordered[,'Protein_number'],ylim=c(8000,9500),las=1)
dev.off()
library(ggplot2) # version 3.3.5
nums <- c(as.numeric(num_tumor_ordered[,'Protein_number']),as.numeric(num_normal_ordered[,'Protein_number']))
nums_df <- data.frame(Numbers=nums,
                      sample_group=c(rep('tumor',103),rep('normal',103)))
pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics/protein_numbers_comp.pdf',height=6,width=5)
ggplot(nums_df,aes(x=sample_group,y=Numbers)) +
  geom_jitter(width=0.15,size=0.7) +
  geom_boxplot(color="black",fill=NA,outlier.shape=NA,width=0.3) +
  scale_x_discrete(limit=c('normal','tumor'),expand=c(0.1,0)) +
  theme(axis.text=element_text(size=10,colour="black"),axis.title=element_blank(),
        panel.background = element_rect(fill = NA,color="black"),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x=element_text(size=7,colour="black",angle=90,vjust=0.5,hjust=1),
        axis.text.y=element_text(size=7,colour="black"),
        strip.text.y = element_text(angle = 0,hjust=0,color="black",size=8),
        strip.background = element_blank(),
        strip.text.x = element_text(color="black",size=7,vjust=0))
dev.off()

### The number distribution of phosphorylation sites in each sample, Fig. 1d
setwd('/home/shengli/projects/AEG_proteomics/results/phosphoproteome')
phosp_num <- read.table('Phosphoproteomics_sample_numbers.txt',header=T,sep='\t')
num_normal <- phosp_num[which(phosp_num[,'Type']=='normal'),]
num_tumor <- phosp_num[which(phosp_num[,'Type']=='tumor'),]
rownames(num_normal) <- num_normal[,'Sample']
rownames(num_tumor) <- num_tumor[,'Sample']
num_normal_ordered <- num_normal[order(num_normal[,'Phosphoproteomics_number'],decreasing = F),]
num_tumor_ordered <- num_tumor[order(num_tumor[,'Phosphoproteomics_number'],decreasing = F),]
pdf('/home/shengli/projects/AEG_proteomics/figures/Phosphoproteomics/Phospho_numbers_tumor.pdf',height=5,width=5)
plot(num_tumor_ordered[,'Phosphoproteomics_number'],ylim=c(5500,10000),las=1)
dev.off()

num_normal <- phosp_num[which(phosp_num[,'Type']=='normal'),]
rownames(num_normal) <- num_normal[,'Sample']
num_normal_ordered <- num_normal[rownames(num_tumor_ordered),]
pdf('/home/shengli/projects/AEG_proteomics/figures/Phosphoproteomics/Protein_numbers_normal.pdf',height=5,width=5)
plot(num_normal_ordered[,'Protein_number'],ylim=c(8000,9500),las=1)
dev.off()

library(ggplot2) # version 3.3.5
nums <- c(as.numeric(num_tumor_ordered[,'Phosphoproteomics_number']),as.numeric(num_normal_ordered[,'Phosphoproteomics_number']))
nums_df <- data.frame(Numbers=nums,
                      sample_group=c(rep('tumor',103),rep('normal',103)))
pdf('/home/shengli/projects/AEG_proteomics/figures/Phosphoproteomics/phosp_numbers_comp.pdf',height=6,width=5)
ggplot(nums_df,aes(x=sample_group,y=Numbers)) +
  geom_jitter(width=0.15,size=0.7) +
  geom_boxplot(color="black",fill=NA,outlier.shape=NA,width=0.3) +
  scale_x_discrete(limit=c('normal','tumor'),expand=c(0.1,0)) +
  theme(axis.text=element_text(size=10,colour="black"),axis.title=element_blank(),
        panel.background = element_rect(fill = NA,color="black"),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x=element_text(size=7,colour="black",angle=90,vjust=0.5,hjust=1),
        axis.text.y=element_text(size=7,colour="black"),
        strip.text.y = element_text(angle = 0,hjust=0,color="black",size=8),
        strip.background = element_blank(),
        strip.text.x = element_text(color="black",size=7,vjust=0))
dev.off()

### The number distribution of expressed genes in each sample, Fig. 1e
setwd('/home/shengli/projects/AEG_proteomics/results/RNAseq')
gene_num <- read.table('gene_sample_numbers.txt',header=T,sep='\t')

num_normal <- gene_num[which(gene_num[,'Type']=='normal'),]
rownames(num_normal) <- num_normal[,'Sample']
num_normal_ordered <- num_normal[order(num_normal[,'Gene_number'],decreasing = F),]
pdf('/home/shengli/projects/AEG_proteomics/figures/RNAseq/Gene_numbers_normal.pdf',height=5,width=5)
barplot(num_normal_ordered[,'Gene_number'],ylim=c(0,30000),las=1)
dev.off()

num_tumor <- gene_num[which(gene_num[,'Type']=='tumor'),]
rownames(num_tumor) <- num_tumor[,'Sample']
num_tumor_ordered <- num_tumor[rownames(num_normal_ordered),]
pdf('/home/shengli/projects/AEG_proteomics/figures/RNAseq/Gene_numbers_tumor.pdf',height=5,width=5)
barplot(num_tumor_ordered[,'Gene_number'],ylim=c(0,30000),las=1)
dev.off()

### Box plots showing the comparisons of TMB between patients with different clinicopathological features, Supplementary Fig. 2
setwd('/home/shengli/projects/AEG_proteomics/data/clinical')
clinical_neat <- read.table('Clinical_infor_neat.txt',header=T,row.names=1,sep='\t')
clinical_neat[which(clinical_neat[,'Age'] < 65),'Age'] <- 1
clinical_neat[which(clinical_neat[,'Age'] >= 65),'Age'] <- 2
clinical_age <- cbind(rownames(clinical_neat),clinical_neat[,'Age'],rep(1,nrow(clinical_neat)))
colnames(clinical_age) <- c('Patient','Age','num')
clinical_gender <- cbind(rownames(clinical_neat),clinical_neat[,'Sex'],rep(1,nrow(clinical_neat)))
colnames(clinical_gender) <- c('Patient','Gender','num')
clinical_smoking <- cbind(rownames(clinical_neat),clinical_neat[,'Smoking'],rep(1,nrow(clinical_neat)))
colnames(clinical_smoking) <- c('Patient','Smoking','num')
clinical_alcohol <- cbind(rownames(clinical_neat),clinical_neat[,'Alcohol'],rep(1,nrow(clinical_neat)))
colnames(clinical_alcohol) <- c('Patient','Alcohol','num')
clinical_siewert <- cbind(rownames(clinical_neat),clinical_neat[,'Siewert.type'],rep(1,nrow(clinical_neat)))
colnames(clinical_siewert) <- c('Patient','Siewert','num')
clinical_stage <- cbind(rownames(clinical_neat),clinical_neat[,'Clinical.stage'],rep(1,nrow(clinical_neat)))
colnames(clinical_stage) <- c('Patient','Stage','num')

setwd('/home/shengli/projects/AEG_proteomics/results/WES')
tumor_tmb <- read.table('AEG_TMB.txt',header=T,row.names=1,sep='\t')
# age
samples_age1 <- clinical_age[which(clinical_age[,'Age']==1),'Patient']
samples_age2 <- clinical_age[which(clinical_age[,'Age']==2),'Patient']
tmb_age1 <- as.numeric(tumor_tmb[samples_age1,'total_perMB'])
tmb_age2 <- as.numeric(tumor_tmb[samples_age2,'total_perMB'])
tmb_age_df <- data.frame(TMB = c(tmb_age1,tmb_age2),
                         Group = c(rep("Group1",length(tmb_age1)),rep("Group2",length(tmb_age2))))
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/TMB_age_boxplot.pdf',width=5,height = 4)
ggplot(tmb_age_df,aes(x=Group,y=TMB)) +
  geom_jitter(width=0.15,size=1.5) +
  geom_boxplot(color="black",fill=NA,outlier.shape=NA,width=0.3) +
  scale_x_discrete(limit=c("Group1","Group2"),expand=c(0.1,0)) +
  theme(axis.text=element_text(size=10,colour="black"),axis.title=element_blank(),
        panel.background = element_rect(fill = NA,color="black"),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x=element_text(size=7,colour="black",angle=90,vjust=0.5,hjust=1),
        axis.text.y=element_text(size=7,colour="black"),
        strip.text.y = element_text(angle = 0,hjust=0,color="black",size=8),
        strip.background = element_blank(),
        strip.text.x = element_text(color="black",size=7,vjust=0))
dev.off()
# sex
samples_sex1 <- clinical_gender[which(clinical_gender[,'Gender']==0),'Patient']
samples_sex2 <- clinical_gender[which(clinical_gender[,'Gender']==1),'Patient']
tmb_sex1 <- as.numeric(tumor_tmb[samples_sex1,'total_perMB'])
tmb_sex2 <- as.numeric(tumor_tmb[samples_sex2,'total_perMB'])
tmb_sex_df <- data.frame(TMB = c(tmb_sex1,tmb_sex2),
                         Group = c(rep("Group1",length(tmb_sex1)),rep("Group2",length(tmb_sex2))))
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/TMB_sex_boxplot.pdf',width=5,height = 4)
ggplot(tmb_sex_df,aes(x=Group,y=TMB)) +
  geom_jitter(width=0.15,size=1.5) +
  geom_boxplot(color="black",fill=NA,outlier.shape=NA,width=0.3) +
  scale_x_discrete(limit=c("Group1","Group2"),expand=c(0.1,0)) +
  theme(axis.text=element_text(size=10,colour="black"),axis.title=element_blank(),
        panel.background = element_rect(fill = NA,color="black"),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x=element_text(size=7,colour="black",angle=90,vjust=0.5,hjust=1),
        axis.text.y=element_text(size=7,colour="black"),
        strip.text.y = element_text(angle = 0,hjust=0,color="black",size=8),
        strip.background = element_blank(),
        strip.text.x = element_text(color="black",size=7,vjust=0))
dev.off()
# smoking
samples_smk1 <- clinical_smoking[which(clinical_smoking[,'Smoking']==0),'Patient']
samples_smk2 <- clinical_smoking[which(clinical_smoking[,'Smoking']==1),'Patient']
tmb_smk1 <- as.numeric(tumor_tmb[samples_smk1,'total_perMB'])
tmb_smk2 <- as.numeric(tumor_tmb[samples_smk2,'total_perMB'])
tmb_smk_df <- data.frame(TMB = c(tmb_smk1,tmb_smk2),
                         Group = c(rep("Group1",length(tmb_smk1)),rep("Group2",length(tmb_smk2))))
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/TMB_smoking_boxplot.pdf',width=5,height = 4)
ggplot(tmb_smk_df,aes(x=Group,y=TMB)) +
  geom_jitter(width=0.15,size=1.5) +
  geom_boxplot(color="black",fill=NA,outlier.shape=NA,width=0.3) +
  scale_x_discrete(limit=c("Group1","Group2"),expand=c(0.1,0)) +
  theme(axis.text=element_text(size=10,colour="black"),axis.title=element_blank(),
        panel.background = element_rect(fill = NA,color="black"),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x=element_text(size=7,colour="black",angle=90,vjust=0.5,hjust=1),
        axis.text.y=element_text(size=7,colour="black"),
        strip.text.y = element_text(angle = 0,hjust=0,color="black",size=8),
        strip.background = element_blank(),
        strip.text.x = element_text(color="black",size=7,vjust=0))
dev.off()
# alcohol
samples_alc1 <- clinical_alcohol[which(clinical_alcohol[,'Alcohol']==0),'Patient']
samples_alc2 <- clinical_alcohol[which(clinical_alcohol[,'Alcohol']==1),'Patient']
tmb_alc1 <- as.numeric(tumor_tmb[samples_alc1,'total_perMB'])
tmb_alc2 <- as.numeric(tumor_tmb[samples_alc2,'total_perMB'])
tmb_alc_df <- data.frame(TMB = c(tmb_alc1,tmb_alc2),
                         Group = c(rep("Group1",length(tmb_alc1)),rep("Group2",length(tmb_alc2))))
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/TMB_alcohol_boxplot.pdf',width=5,height = 4)
ggplot(tmb_alc_df,aes(x=Group,y=TMB)) +
  geom_jitter(width=0.15,size=1.5) +
  geom_boxplot(color="black",fill=NA,outlier.shape=NA,width=0.3) +
  scale_x_discrete(limit=c("Group1","Group2"),expand=c(0.1,0)) +
  theme(axis.text=element_text(size=10,colour="black"),axis.title=element_blank(),
        panel.background = element_rect(fill = NA,color="black"),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x=element_text(size=7,colour="black",angle=90,vjust=0.5,hjust=1),
        axis.text.y=element_text(size=7,colour="black"),
        strip.text.y = element_text(angle = 0,hjust=0,color="black",size=8),
        strip.background = element_blank(),
        strip.text.x = element_text(color="black",size=7,vjust=0))
dev.off()
# Siewert type
samples_swt1 <- clinical_siewert[which(clinical_siewert[,'Siewert']==1),'Patient']
samples_swt2 <- clinical_siewert[which(clinical_siewert[,'Siewert']==2),'Patient']
samples_swt3 <- clinical_siewert[which(clinical_siewert[,'Siewert']==3),'Patient']
tmb_swt1 <- as.numeric(tumor_tmb[samples_swt1,'total_perMB'])
tmb_swt2 <- as.numeric(tumor_tmb[samples_swt2,'total_perMB'])
tmb_swt3 <- as.numeric(tumor_tmb[samples_swt3,'total_perMB'])

tmb_swt_df <- data.frame(TMB = c(tmb_swt1,tmb_swt2,tmb_swt3),
                         Group = c(rep("Group1",length(tmb_swt1)),rep("Group2",length(tmb_swt2)),rep("Group3",length(tmb_swt3))))

pdf('/home/shengli/projects/AEG_proteomics/figures/WES/TMB_siewert_boxplot.pdf',width=5,height = 4)
ggplot(tmb_swt_df,aes(x=Group,y=TMB)) +
  geom_jitter(width=0.15,size=1.5) +
  geom_boxplot(color="black",fill=NA,outlier.shape=NA,width=0.3) +
  scale_x_discrete(limit=c("Group1","Group2","Group3"),expand=c(0.1,0)) +
  theme(axis.text=element_text(size=10,colour="black"),axis.title=element_blank(),
        panel.background = element_rect(fill = NA,color="black"),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x=element_text(size=7,colour="black",angle=90,vjust=0.5,hjust=1),
        axis.text.y=element_text(size=7,colour="black"),
        strip.text.y = element_text(angle = 0,hjust=0,color="black",size=8),
        strip.background = element_blank(),
        strip.text.x = element_text(color="black",size=7,vjust=0))
dev.off()
# Stage
samples_stg1 <- clinical_stage[which(clinical_stage[,'Stage']==1),'Patient']
samples_stg2 <- clinical_stage[which(clinical_stage[,'Stage']==2),'Patient']
samples_stg3 <- clinical_stage[which(clinical_stage[,'Stage']==3),'Patient']
samples_stg4 <- clinical_stage[which(clinical_stage[,'Stage']==4),'Patient']
tmb_stg1 <- as.numeric(tumor_tmb[samples_stg1,'total_perMB'])
tmb_stg2 <- as.numeric(tumor_tmb[samples_stg2,'total_perMB'])
tmb_stg3 <- as.numeric(tumor_tmb[samples_stg3,'total_perMB'])
tmb_stg4 <- as.numeric(tumor_tmb[samples_stg4,'total_perMB'])
tmb_stg_df <- data.frame(TMB = c(tmb_stg1,tmb_stg2,tmb_stg3,tmb_stg4),
                         Group = c(rep("Group1",length(tmb_stg1)),rep("Group2",length(tmb_stg2)),rep("Group3",length(tmb_stg3)),rep("Group4",length(tmb_stg4))))
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/TMB_stage_boxplot.pdf',width=5,height = 4)
ggplot(tmb_stg_df,aes(x=Group,y=TMB)) +
  geom_jitter(width=0.15,size=1.5) +
  geom_boxplot(color="black",fill=NA,outlier.shape=NA,width=0.3) +
  scale_x_discrete(limit=c("Group1","Group2","Group3","Group4"),expand=c(0.1,0)) +
  theme(axis.text=element_text(size=10,colour="black"),axis.title=element_blank(),
        panel.background = element_rect(fill = NA,color="black"),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x=element_text(size=7,colour="black",angle=90,vjust=0.5,hjust=1),
        axis.text.y=element_text(size=7,colour="black"),
        strip.text.y = element_text(angle = 0,hjust=0,color="black",size=8),
        strip.background = element_blank(),
        strip.text.x = element_text(color="black",size=7,vjust=0))
dev.off()

### Distribution of relative abundance of detected proteins, Supplementary Fig. 3a
rm(list=ls())
library(ggplot2)
setwd('/home/shengli/projects/AEG_proteomics/data/proteome')
prot_ibaq <- read.table('Proteomics_iBAQ103_log2quantile_normlization_df.txt',header=T)
sample_cancer <- as.character(unique(prot_ibaq[which(prot_ibaq[,'Group']=='Cancer'),'Sample']))
sample_normal <- as.character(unique(prot_ibaq[which(prot_ibaq[,'Group']=='Normal'),'Sample']))
samples <- c(sample_cancer,sample_normal)

pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics_iBAQ_norm_distribution.pdf',width=20,height=5)
ggplot(prot_ibaq,aes(x=Sample,y=Abundance,fill=Group)) +
  geom_boxplot()
dev.off()

### Boxplot showing the comparison of the numbers of expressed genes between parired AEG tumor and NAT samples, Supplementary Fig. 4c
library(ggplot2) # version 3.3.5
nums <- c(as.numeric(num_tumor_ordered[,'Gene_number']),as.numeric(num_normal_ordered[,'Gene_number']))
nums_df <- data.frame(Numbers=nums,
                      sample_group=c(rep('tumor',85),rep('normal',85)))
pdf('/home/shengli/projects/AEG_proteomics/figures/RNAseq/gene_numbers_comp.pdf',height=6,width=5)
ggplot(nums_df,aes(x=sample_group,y=Numbers)) +
  geom_jitter(width=0.15,size=0.7) +
  geom_boxplot(color="black",fill=NA,outlier.shape=NA,width=0.3) +
  scale_x_discrete(limit=c('normal','tumor'),expand=c(0.1,0)) +
  theme(axis.text=element_text(size=10,colour="black"),axis.title=element_blank(),
        panel.background = element_rect(fill = NA,color="black"),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x=element_text(size=7,colour="black",angle=90,vjust=0.5,hjust=1),
        axis.text.y=element_text(size=7,colour="black"),
        strip.text.y = element_text(angle = 0,hjust=0,color="black",size=8),
        strip.background = element_blank(),
        strip.text.x = element_text(color="black",size=7,vjust=0))
dev.off()
wilcox.test(as.numeric(nums_df[which(nums_df[,'sample_group']=='normal'),'Numbers']),as.numeric(nums_df[which(nums_df[,'sample_group']=='tumor'),'Numbers']))


