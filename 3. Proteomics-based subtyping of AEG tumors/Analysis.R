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

## differential protein expression among different subtypes
library(limma)
setwd('/home/shengli/projects/AEG_proteomics/data')
groups <- read.table('survival_data_v3.txt',header=T,sep='\t')
samples_g1 <- groups[which(groups[,'Group']=='Group1'),'Patient']
samples_g2 <- groups[which(groups[,'Group']=='Group2'),'Patient']
samples_g3 <- groups[which(groups[,'Group']=='Group3'),'Patient']
setwd('/home/shengli/projects/AEG_proteomics/results')
prot_mx <- read.table('Proteomics_subgroups_perc25.txt',header=T,row.names=1,sep='\t')
g1_g2_mx <- prot_mx[,c(samples_g1,samples_g2)]
g1_g3_mx <- prot_mx[,c(samples_g1,samples_g3)]
g2_g3_mx <- prot_mx[,c(samples_g2,samples_g3)]

design_g1g2 <- cbind(Intercept=1,
                     Group=c(rep(1,length(samples_g1)),rep(0,length(samples_g2))))
fit <- lmFit(g1_g2_mx,design_g1g2)
fit <- eBayes(fit,trend = TRUE)
deg_g1g2 <- topTable(fit,coef=2,number=Inf)
setwd('/home/shengli/projects/AEG_proteomics/results')
write.table(deg_g1g2,file='G1_G2_DEP_all.txt',quote=F,sep='\t',append=F)

design_g1g3 <- cbind(Intercept=1,
                     Group=c(rep(1,length(samples_g1)),rep(0,length(samples_g3))))
fit <- lmFit(g1_g3_mx,design_g1g3)
fit <- eBayes(fit,trend = TRUE)
deg_g1g3 <- topTable(fit,coef=2,number=Inf)
setwd('/home/shengli/projects/AEG_proteomics/results')
write.table(deg_g1g3,file='G1_G3_DEP_all.txt',quote=F,sep='\t',append=F)

design_g2g3 <- cbind(Intercept=1,
                     Group=c(rep(1,length(samples_g2)),rep(0,length(samples_g3))))
fit <- lmFit(g2_g3_mx,design_g2g3)
fit <- eBayes(fit,trend = TRUE)
deg_g2g3 <- topTable(fit,coef=2,number=Inf)
setwd('/home/shengli/projects/AEG_proteomics/results')
write.table(deg_g2g3,file='G2_G3_DEP_all.txt',quote=F,sep='\t',append=F)

## differential proteins between subtype tumor and paired normal samples
library(limma)
setwd('/home/shengli/projects/AEG_proteomics/data')
groups <- read.table('./clinical/survival_data_v3.txt',header=T,sep='\t')
samples_g1 <- groups[which(groups[,'Group']=='Group1'),'Patient']
samples_g2 <- groups[which(groups[,'Group']=='Group2'),'Patient']
samples_g3 <- groups[which(groups[,'Group']=='Group3'),'Patient']
mx_g1 <- read.table('./proteome/G1_tumor_nat_perc25.txt',header=T,row.names=1,sep='\t')
design_g1 <- cbind(Intercept=1,
                   Group=c(rep(1,length(samples_g1)),rep(0,length(samples_g1))))
fit <- lmFit(mx_g1,design_g1)
fit <- eBayes(fit,trend = TRUE)
deg_g1 <- topTable(fit,coef=2,number=Inf)
setwd('/home/shengli/projects/AEG_proteomics/results')
write.table(deg_g1,file='G1_NAT_DEP_all.txt',quote=F,sep='\t',append=F)

setwd('/home/shengli/projects/AEG_proteomics/data')
mx_g2 <- read.table('G2_tumor_nat_perc25.txt',header=T,row.names=1,sep='\t')
design_g2 <- cbind(Intercept=1,
                   Group=c(rep(1,length(samples_g2)),rep(0,length(samples_g2))))
fit <- lmFit(mx_g2,design_g2)
fit <- eBayes(fit,trend = TRUE)
deg_g2 <- topTable(fit,coef=2,number=Inf)
setwd('/home/shengli/projects/AEG_proteomics/results')
write.table(deg_g2,file='G2_NAT_DEP_all.txt',quote=F,sep='\t',append=F)

setwd('/home/shengli/projects/AEG_proteomics/data')
mx_g3 <- read.table('G3_tumor_nat_perc25.txt',header=T,row.names=1,sep='\t')
design_g3 <- cbind(Intercept=1,
                   Group=c(rep(1,length(samples_g3)),rep(0,length(samples_g3))))
fit <- lmFit(mx_g3,design_g3)
fit <- eBayes(fit,trend = TRUE)
deg_g3 <- topTable(fit,coef=2,number=Inf)
setwd('/home/shengli/projects/AEG_proteomics/results')
write.table(deg_g3,file='G3_NAT_DEP_all.txt',quote=F,sep='\t',append=F)

## identify signature proteins for each subtype
setwd('/home/shengli/projects/AEG_proteomics/data')
groups <- read.table('survival_data_v3.txt',header=T,sep='\t')
samples_g1 <- groups[which(groups[,'Group']=='Group1'),'Patient']
samples_g2 <- groups[which(groups[,'Group']=='Group2'),'Patient']
samples_g3 <- groups[which(groups[,'Group']=='Group3'),'Patient']
mx_all <- read.table('Proteomics_iBAQ_log2quantile_normlization_impute.txt',header=T,row.names=1,sep='\t')
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
dep_g1_nat <- read.table('G1_NAT_DEP_all.txt',header=T,row.names=1,sep='\t')
dep_g2_nat <- read.table('G2_NAT_DEP_all.txt',header=T,row.names=1,sep='\t')
dep_g3_nat <- read.table('G3_NAT_DEP_all.txt',header=T,row.names=1,sep='\t')
dep_g1_g2 <- read.table('G1_G2_DEP_all.txt',header=T,row.names=1,sep='\t')
dep_g1_g3 <- read.table('G1_G3_DEP_all.txt',header=T,row.names=1,sep='\t')
dep_g2_g3 <- read.table('G2_G3_DEP_all.txt',header=T,row.names=1,sep='\t')

# G1 subtype
dep_g1_up <- rownames(dep_g1_nat)[which(dep_g1_nat[,'adj.P.Val'] < 0.05 & dep_g1_nat[,'logFC'] > 1)]
dep_g1g2_up <- rownames(dep_g1_g2)[which(dep_g1_g2[,'adj.P.Val'] < 0.05 & dep_g1_g2[,'logFC'] < -1)]
dep_g1g3_up <- rownames(dep_g1_g3)[which(dep_g1_g3[,'adj.P.Val'] < 0.05 & dep_g1_g3[,'logFC'] < -1)]
sigp_g1 <- intersect(intersect(dep_g1_up,dep_g1g2_up),dep_g1g3_up)
# G2 subtype
dep_g2_up <- rownames(dep_g2_nat)[which(dep_g2_nat[,'adj.P.Val'] < 0.05 & dep_g2_nat[,'logFC'] > 1)]
dep_g1g2_up <- rownames(dep_g1_g2)[which(dep_g1_g2[,'adj.P.Val'] < 0.05 & dep_g1_g2[,'logFC'] > 1)]
dep_g2g3_up <- rownames(dep_g2_g3)[which(dep_g2_g3[,'adj.P.Val'] < 0.05 & dep_g2_g3[,'logFC'] < -1)]
sigp_g2 <- intersect(intersect(dep_g2_up,dep_g1g2_up),dep_g2g3_up)
# G3 subtype
dep_g3_up <- rownames(dep_g3_nat)[which(dep_g3_nat[,'adj.P.Val'] < 0.05 & dep_g3_nat[,'logFC'] > 1)]
dep_g1g3_up <- rownames(dep_g1_g3)[which(dep_g1_g3[,'adj.P.Val'] < 0.05 & dep_g1_g3[,'logFC'] > 1)]
dep_g2g3_up <- rownames(dep_g2_g3)[which(dep_g2_g3[,'adj.P.Val'] < 0.05 & dep_g2_g3[,'logFC'] > 1)]
sigp_g3 <- intersect(intersect(dep_g3_up,dep_g1g3_up),dep_g2g3_up)

sigp_protein <- c(sigp_g1,sigp_g2,sigp_g3)
sigp_subtypes <- c(rep('Group1',length(sigp_g1)),rep('Group2',length(sigp_g2)),rep('Group3',length(sigp_g3)))
sigp <- cbind(sigp_protein,sigp_subtypes)
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
write.table(sigp,file='Subtype_signature_proteins.txt',sep='\t',quote=F,append=F)

# mutation comparison between subtypes
setwd('/home/shengli/projects/AEG_proteomics/results/WES')
num_g1 <- 40
num_g2 <- 23
num_g3 <- 40
muts_gp <- read.table('AEG_groups_mutation.txt',header=T,sep='\t',row.names=1)
results_g1_pvalues <- c()
results_g1_ors <- c()
results_g2_pvalues <- c()
results_g2_ors <- c()
results_g3_pvalues <- c()
results_g3_ors <- c()
for (n in 1:nrow(muts_gp)) {
  g1_mut_sub <- as.numeric(muts_gp[n,1])
  g1_mut_nonsub <- as.numeric(muts_gp[n,2]) + as.numeric(muts_gp[n,3])
  g1_nonmut_sub <- num_g1 - as.numeric(muts_gp[n,1])
  g1_nonmut_nonsub <- num_g2 + num_g3 - as.numeric(muts_gp[n,2]) - as.numeric(muts_gp[n,3])
  g1_mat <- matrix(c(g1_mut_sub,g1_mut_nonsub,g1_nonmut_sub,g1_nonmut_nonsub),nrow=2)
  g1_test <- fisher.test(g1_mat)
  results_g1_ors <- c(results_g1_ors,as.numeric(g1_test$estimate))
  results_g1_pvalues <- c(results_g1_pvalues,as.numeric(g1_test$p.value))
  g2_mut_sub <- as.numeric(muts_gp[n,2])
  g2_mut_nonsub <- as.numeric(muts_gp[n,1]) + as.numeric(muts_gp[n,3])
  g2_nonmut_sub <- num_g2 - as.numeric(muts_gp[n,2])
  g2_nonmut_nonsub <- num_g1 + num_g3 - as.numeric(muts_gp[n,1]) - as.numeric(muts_gp[n,3])
  g2_mat <- matrix(c(g2_mut_sub,g2_mut_nonsub,g2_nonmut_sub,g2_nonmut_nonsub),nrow=2)
  g2_test <- fisher.test(g2_mat)
  results_g2_ors <- c(results_g2_ors,as.numeric(g2_test$estimate))
  results_g2_pvalues <- c(results_g2_pvalues,as.numeric(g2_test$p.value))
  g3_mut_sub <- as.numeric(muts_gp[n,3])
  g3_mut_nonsub <- as.numeric(muts_gp[n,1]) + as.numeric(muts_gp[n,2])
  g3_nonmut_sub <- num_g3 - as.numeric(muts_gp[n,3])
  g3_nonmut_nonsub <- num_g1 + num_g2 - as.numeric(muts_gp[n,1]) - as.numeric(muts_gp[n,2])
  g3_mat <- matrix(c(g3_mut_sub,g3_mut_nonsub,g3_nonmut_sub,g3_nonmut_nonsub),nrow=2)
  g3_test <- fisher.test(g3_mat)
  results_g3_ors <- c(results_g3_ors,as.numeric(g3_test$estimate))
  results_g3_pvalues <- c(results_g3_pvalues,as.numeric(g3_test$p.value))
}

results <- cbind(as.character(rownames(muts_gp)),results_g1_ors,results_g1_pvalues,results_g2_ors,results_g2_pvalues,results_g3_ors,results_g3_pvalues)
colnames(results) <- c('Gene','Group1_OR','Group1_Pvalue','Group2_OR','Group2_Pvalue','Group3_OR','Group3_Pvalue')
setwd('/home/shengli/projects/AEG_proteomics/results/WES')
write.table(results,file='subgroups_mut_fishertest.txt',sep='\t',quote=F,row.names=F,append=F)

### transcriptional activities of hallmark comparisons between subtypes
rm(list=ls())
library(GSVA)
hallmarks <- read.table('/home/public/public_data/MSigDB/processed/h.all.v6.2.symbols.txt',header=T,sep='\t',row.names=1)
dir_expr <- '/home/shengli/projects/AEG_proteomics/results/RNAseq'
file_aeg <- paste(dir_expr,'/','AEG_gene_tmp_mx_ordered.txt',sep='')
expr_aeg <- read.table(file_aeg,header=T,row.names=1,sep='\t')
rownames(expr_aeg) <- toupper(rownames(expr_aeg))
expr_aeg <- as.matrix(expr_aeg)
expr_aeg <- log2(1+expr_aeg)
es_aeg <- matrix(ncol=ncol(expr_aeg),nrow=nrow(hallmarks))
for (n in 1:nrow(hallmarks)) {
  hallmark_genes <- unlist(strsplit(as.character(hallmarks[n,1]),'[|]'))
  gene_list <- list(hallmark=hallmark_genes)
  es_hk <- gsva(expr_aeg,gene_list,mx.diff=FALSE,method='ssgsea',verbose=FALSE,parallel.sz=1)
  es_aeg[n,] <- es_hk
}
rownames(es_aeg) <- rownames(hallmarks)
colnames(es_aeg) <- colnames(expr_aeg)
setwd('/home/shengli/projects/AEG_proteomics/results/RNAseq')
write.table(es_aeg,file='AEG_hallmark_scores.txt',quote=F,sep='\t')

# difference of gene-level hallmarks scores between subtypes
results_pvalues_s1s2 <- c()
results_foldchanges_s1s2 <- c()
results_pvalues_s1s3 <- c()
results_foldchanges_s1s3 <- c()
results_pvalues_s2s3 <- c()
results_foldchanges_s2s3 <- c()
for (n in 1:nrow(hallmark_scores)) {
  hk_scores_s1 <- as.numeric(hallmark_scores[n,samples_t1_p])
  hk_scores_s2 <- as.numeric(hallmark_scores[n,samples_t2_p])
  hk_scores_s3 <- as.numeric(hallmark_scores[n,samples_t3_p])
  ave_hk_s1 <- mean(hk_scores_s1)
  ave_hk_s2 <- mean(hk_scores_s2)
  ave_hk_s3 <- mean(hk_scores_s3)
  ave_hk_fc_s1s2 <- round(ave_hk_s2/ave_hk_s1,3)
  results_foldchanges_s1s2 <- c(results_foldchanges_s1s2,ave_hk_fc_s1s2)
  test_hk_s1s2 <- t.test(hk_scores_s1,hk_scores_s2)
  results_pvalues_s1s2 <- c(results_pvalues_s1s2,test_hk_s1s2$p.value)
  ave_hk_fc_s1s3 <- round(ave_hk_s3/ave_hk_s1,3)
  results_foldchanges_s1s3 <- c(results_foldchanges_s1s3,ave_hk_fc_s1s3)
  test_hk_s1s3 <- t.test(hk_scores_s1,hk_scores_s3)
  results_pvalues_s1s3 <- c(results_pvalues_s1s3,test_hk_s1s3$p.value)
  ave_hk_fc_s2s3 <- round(ave_hk_s3/ave_hk_s2,3)
  results_foldchanges_s2s3 <- c(results_foldchanges_s2s3,ave_hk_fc_s2s3)
  test_hk_s2s3 <- t.test(hk_scores_s2,hk_scores_s3)
  results_pvalues_s2s3 <- c(results_pvalues_s2s3,test_hk_s2s3$p.value)
}
results_fdr_s1s2 <- p.adjust(results_pvalues_s1s2,method='fdr')
results_fdr_s1s3 <- p.adjust(results_pvalues_s1s3,method='fdr')
results_fdr_s2s3 <- p.adjust(results_pvalues_s2s3,method='fdr')
results <- cbind(results_foldchanges_s1s2,results_pvalues_s1s2,results_fdr_s1s2,results_foldchanges_s1s3,results_pvalues_s1s3,results_fdr_s1s3,
                 results_foldchanges_s2s3,results_pvalues_s2s3,results_fdr_s2s3)
rownames(results) <- rownames(hallmark_scores)
setwd('/home/shengli/projects/AEG_proteomics/results/RNAseq')
write.table(results,file='AEG_subtypes_hallmarks_diff.txt',sep='\t',quote=F,append=F)

