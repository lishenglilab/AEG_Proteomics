## plot individual survival curves of proteins, Supplementary Fig. 10
library(survival)
library(ggplot2)
library(survminer)
setwd('/home/shengli/projects/AEG_proteomics/data/clinical')
surv_data <- read.table('survival_data_v3.txt',header=T,sep='\t')
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
prot_profile <- read.table('Proteomics_iBAQ103_log2quantile_normlization_impute_perc25.txt',header=T,row.names=1,sep='\t')
prot_profile_tumor <- prot_profile['P42771',1:103]
mid_prot <- median(as.numeric(prot_profile_tumor))
pt_high <- colnames(prot_profile_tumor)[which(prot_profile_tumor[1,] > mid_prot)]
pt_low <- colnames(prot_profile_tumor)[which(prot_profile_tumor[1,] <= mid_prot)]
surv_data[which(surv_data[,'Patient'] %in% pt_high),'Group'] <- 'High'
surv_data[which(surv_data[,'Patient'] %in% pt_low),'Group'] <- 'Low'
fit <- survfit(Surv(Time,Status)~Group, data=surv_data)
# summary(cox_model) 0. (95% CI, 0.35-0.40)
pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics/P35869_survival.pdf',width=5.33,height=4.82)
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = surv_data,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = FALSE,         # show confidence intervals for 
  xlim = c(0,100),         # present narrower X axis, but not affect
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 20,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
)
dev.off()

# plot significant mutations between subtypes, Fig. 3c
rm(list=ls())
library(ggplot2)
library(ggrepel)
setwd('/home/shengli/projects/AEG_proteomics/results/WES')
enrich <- read.table('subgroups_mut_sig.txt',header=T,sep='\t',row.names=1)
enrich$lgp <- -log10(enrich[,'Pvalue'])
enrich[which(enrich[,'OR'] == 0),'OR'] <- 1
enrich$logor <- log2(enrich[,'OR'])
enrich$gene <- row.names(enrich)
enrich_none <- enrich[which(enrich[,'Group'] == 'Group0'),]
enrich_g1 <- enrich[which(enrich[,'Group'] == 'Group1'),]
enrich_g2 <- enrich[which(enrich[,'Group'] == 'Group2'),]
enrich_g3 <- enrich[which(enrich[,'Group'] == 'Group3'),]
g1_gene <- as.character(row.names(enrich_g1))
g2_gene <- as.character(row.names(enrich_g2))
g3_gene <- as.character(row.names(enrich_g3))
freq <- read.table('AEG_groups_mutation.txt',header=T,sep='\t')
freq_g1 <- freq[which(freq[,'Gene'] %in% g1_gene),]
freq_g2 <- freq[which(freq[,'Gene'] %in% g2_gene),]
freq_g3 <- freq[which(freq[,'Gene'] %in% g3_gene),]
enrich_g1$freq <- freq_g1[,'Group1']
enrich_g2$freq <- freq_g2[,'Group2']
enrich_g3$freq <- freq_g3[,'Group3']
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/AEG_mutation_enrich_groups.pdf',width=6.87,height=5.92)
ggplot(enrich,aes(x=logor,y=lgp,label=gene)) +
  geom_point(data=enrich_none,aes(x=logor,y=lgp),color='lightgray') +
  geom_point(data=enrich_g1,aes(x=logor,y=lgp,size=freq),color='red') +
  geom_point(data=enrich_g2,aes(x=logor,y=lgp,size=freq),color='blue') +
  geom_point(data=enrich_g3,aes(x=logor,y=lgp,size=freq),color='green') +
  geom_label_repel(data = subset(enrich, lgp > 2 & abs(logor) > 2),
                   size = 3,
                   box.padding = 0.5,
                   point.padding = 0.5,
                   segment.size = 0.2,
                   segment.color = "grey50") +
  theme(panel.background=element_rect(colour="black",fill="white"))
dev.off()

# plot individual genes that have significant mutations between subtypes
setwd('/home/shengli/projects/AEG_proteomics/results/WES')
num_g1 <- 40
num_g2 <- 23
num_g3 <- 40
muts_gp <- read.table('AEG_groups_mutation.txt',header=T,sep='\t',row.names=1)
gene <- 'WIZ'
mut_freq <- muts_gp[gene,]
mut_rate <- c(0/40,0/23,5/40)
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/WIZ_subtype_mutation_freq.pdf',height=5,width=3)
barplot(mut_rate,las=2,ylim=c(0,0.15))
dev.off()


