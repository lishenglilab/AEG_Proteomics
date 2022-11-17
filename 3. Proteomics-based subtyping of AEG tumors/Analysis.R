### survival analysis between AEG subtypes
library(survival)
library(ggplot2)
library(survminer)
setwd('/home/shengli/projects/0.Collaborations/4.Qin.Jiangjiang/data/clinical')
surv_data <- read.table('survival_data_v3.txt',header=T,sep='\t')
surv_12 <- surv_data[which(surv_data[,'Group'] != 'Group3'),]
fit <- survfit(Surv(Time,Status)~Group, data=surv_12)
cox_model <- coxph(Surv(Time,Status)~Group, data=surv_12)
# summary(cox_model) 0.71 (95% CI, 0.35-1.43)
pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics_groups12_survival.pdf',width=5.33,height=4.82) # Fig. 3b
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = surv_12,             # data used to fit survival curves.
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
