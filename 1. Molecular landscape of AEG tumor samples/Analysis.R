### Calculate and compare TMB

library(maftools)
library(ggplot2)
setwd('/home/shengli/projects/0.Collaborations/4.Qin.Jiangjiang/data/WES')
aeg <- read.maf(maf='AEG_mutation_neat_filter_non151058_cancer_gene.maf')

aeg_tmb <- tmb(aeg,captureSize = 1,logScale = F)
setwd('/home/shengli/projects/0.Collaborations/4.Qin.Jiangjiang/results/WES')
write.table(aeg_tmb, file='AEG_TMB.txt',quote=F,sep='\t',row.names=F)

tumor_tmb <- read.table('AEG_TMB.txt',header=T,row.names=1,sep='\t')

# age
samples_age1 <- clinical_age[which(clinical_age[,'Age']==1),'Patient']
samples_age2 <- clinical_age[which(clinical_age[,'Age']==2),'Patient']
tmb_age1 <- as.numeric(tumor_tmb[samples_age1,'total_perMB'])
tmb_age2 <- as.numeric(tumor_tmb[samples_age2,'total_perMB'])
comp_age <- wilcox.test(tmb_age1,tmb_age2) # p-value = 0.04463
# Sex
samples_sex1 <- clinical_gender[which(clinical_gender[,'Gender']==0),'Patient']
samples_sex2 <- clinical_gender[which(clinical_gender[,'Gender']==1),'Patient']
tmb_sex1 <- as.numeric(tumor_tmb[samples_sex1,'total_perMB'])
tmb_sex2 <- as.numeric(tumor_tmb[samples_sex2,'total_perMB'])
comp_sex <- wilcox.test(tmb_sex1,tmb_sex2) # p-value = 0.05509
# Smoking
samples_smk1 <- clinical_smoking[which(clinical_smoking[,'Smoking']==0),'Patient']
samples_smk2 <- clinical_smoking[which(clinical_smoking[,'Smoking']==1),'Patient']
tmb_smk1 <- as.numeric(tumor_tmb[samples_smk1,'total_perMB'])
tmb_smk2 <- as.numeric(tumor_tmb[samples_smk2,'total_perMB'])
comp_smk <- wilcox.test(tmb_smk1,tmb_smk2) # p-value = 0.9757
# Alcohol
samples_alc1 <- clinical_alcohol[which(clinical_alcohol[,'Alcohol']==0),'Patient']
samples_alc2 <- clinical_alcohol[which(clinical_alcohol[,'Alcohol']==1),'Patient']
tmb_alc1 <- as.numeric(tumor_tmb[samples_alc1,'total_perMB'])
tmb_alc2 <- as.numeric(tumor_tmb[samples_alc2,'total_perMB'])
comp_alc <- wilcox.test(tmb_alc1,tmb_alc2) # p-value = 0.1379
# Siewert type
samples_swt1 <- clinical_siewert[which(clinical_siewert[,'Siewert']==1),'Patient']
samples_swt2 <- clinical_siewert[which(clinical_siewert[,'Siewert']==2),'Patient']
samples_swt3 <- clinical_siewert[which(clinical_siewert[,'Siewert']==3),'Patient']
tmb_swt1 <- as.numeric(tumor_tmb[samples_swt1,'total_perMB'])
tmb_swt2 <- as.numeric(tumor_tmb[samples_swt2,'total_perMB'])
tmb_swt3 <- as.numeric(tumor_tmb[samples_swt3,'total_perMB'])

tmb_swt_df <- data.frame(TMB = c(tmb_swt1,tmb_swt2,tmb_swt3),
                         Group = c(rep("Group1",length(tmb_swt1)),rep("Group2",length(tmb_swt2)),rep("Group3",length(tmb_swt3))))
comp_swt <- kruskal.test(TMB~Group,data=tmb_swt_df) # p-value = 0.7241
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
comp_stg <- kruskal.test(TMB~Group,data=tmb_stg_df) # p-value = 0.6934

