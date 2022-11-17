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


