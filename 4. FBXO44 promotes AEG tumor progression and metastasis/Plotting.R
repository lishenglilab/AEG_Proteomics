### boxplot for FBXO44 protein difference, Fig. 4a
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
fbxo44_mx <- read.table('Q9H4M3_FBXO44_groups_df.txt',header=T,sep='\t')
fbxo44_mx_min <- fbxo44_mx[which(fbxo44_mx[,'Abundance'] > 0),]

pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics/Q9H4M3_FBXO44_groups_boxplot.pdf',width=5,height = 4)
ggplot(fbxo44_mx_min,aes(x=Group,y=Abundance)) +
  geom_jitter(width=0.15,size=1.5) +
  geom_boxplot(color="black",fill=NA,outlier.shape=NA,width=0.3) +
  scale_x_discrete(limit=c("Group1_normal","Group1_tumor","Group2_normal","Group2_tumor","Group3_normal","Group3_tumor"),expand=c(0.1,0)) +
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

## expression difference of FBXO44 in TCGA pan-cancer cohorts, Supplementary Fig. 11a
setwd('/home/shengli/projects/AEG_proteomics/data/public_data')
expr <- read.table('FBXO44_tcga_paired_expression.txt',header=T,sep='\t')
groups <- c('Normal_BLCA','Tumor_BLCA','Normal_BRCA','Tumor_BRCA','Normal_CESC','Tumor_CESC','Normal_CHOL','Tumor_CHOL','Normal_COAD','Tumor_COAD',
            'Normal_ESCA','Tumor_ESCA','Normal_HNSC','Tumor_HNSC','Normal_KICH','Tumor_KICH','Normal_KIRC','Tumor_KIRC','Normal_KIRP','Tumor_KIRP',
            'Normal_LIHC','Tumor_LIHC','Normal_LUAD','Tumor_LUAD','Normal_LUSC','Tumor_LUSC','Normal_PRAD','Tumor_PRAD','Normal_READ','Tumor_READ',
            'Normal_STAD','Tumor_STAD','Normal_THCA','Tumor_THCA','Normal_UCEC','Tumor_UCEC')
pdf('/home/shengli/projects/AEG_proteomics/figures/FBXO44_tcga_paired_expression_boxplot.pdf',width=10,height = 5.5)
ggplot(expr,aes(x=Group,y=Expression)) +
  geom_jitter(width=0.15,size=1.5) +
  geom_boxplot(color="black",fill=NA,outlier.shape=NA,width=0.3) +
  scale_x_discrete(limit=groups,expand=c(0.1,0)) +
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

### survival analysis of FBXO44, Fig. 4c
library(survival)
library(ggplot2)
library(survminer)
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
data_sigp <- read.table('sigp_for_coxp_v2.txt',header=T,sep='\t')

protein <- 'Q9H4M3'
abundances <- as.numeric(data_sigp[,protein])
median_abund <- median(abundances)
data_sigp$Group <- rep('High',nrow(data_sigp))
data_sigp[which(data_sigp[,protein] <= median_abund),'Group'] <- 'Low'
fit <- survfit(Surv(Time,Status)~Group, data=data_sigp)
cox_model <- coxph(Surv(Time,Status)~Group+Age+Sex+Smoking+Alcohol+Siewert_type+Clinical_stage, data=data_sigp)
summary(cox_model)
# summary(cox_model) HR=0.48 (95% CI, 0.26-0.88)
pdf('/home/shengli/projects/AEG_proteomics/figures/Survival_curve_Q9H4M3.pdf',height=5,width=5)
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = data_sigp,             # data used to fit survival curves.
  risk.table = FALSE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = FALSE,         # show confidence intervals for 
  xlim = c(0,60),         # present narrower X axis, but not affect
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 10,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
)
dev.off()

###

