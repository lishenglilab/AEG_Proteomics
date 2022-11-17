### normalize the protein abundance iBAQ
library(limma)
setwd('/home/shengli/projects/AEG_proteomics/data/proteome')
prot_ibaq <- read.table('Proteomics_iBAQ103.txt',header=T,row.names=1)
prot_ibaq <- as.matrix(prot_ibaq)
quant_norm_ibaq <- normalizeQuantiles(prot_ibaq)
write.table(quant_norm_ibaq,file='Proteomics_iBAQ103_quantile_normlization.txt',quote=F,sep='\t')
log2_norm_ibaq <- log2(quant_norm_ibaq)
write.table(log2_norm_ibaq,file='Proteomics_iBAQ103_log2quantile_normlization.txt',quote=F,sep='\t')

### differentially expressed proteins analysis
library(limma) # version 3.46.0
setwd('/home/shengli/projects/AEG_proteomics/data/proteome')
prot_mx <- read.table('Proteomics_iBAQ103_log2quantile_normlization_impute_perc25.txt',header=T,row.names=1,sep='\t')
design_nt <- cbind(Intercept=1,
                   Group=c(rep(1,103),rep(0,103)))
fit <- lmFit(prot_mx,design_nt)
fit <- eBayes(fit,trend=TRUE)
dep_nt <- topTable(fit,coef=2,number=Inf)
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
write.table(dep_nt,file="Normal_Tumor_DEP_all.txt",quote=F,sep='\t',append=F)

### enrichment analysis for up-regulated proteins
rm(list = ls())
options(stringsAsFactors = FALSE) # prohibit shift form character to factor
library(tidyverse)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
dep_df <- read.table('Normal_Tumor_DEP_sig_genes.txt',header=T,row.names=1,sep='\t')
up_df <- dep_df[which(dep_df[,'logFC']>0),]
degID <- bitr(up_df$Gene, fromType = "SYMBOL", toType = c( "ENTREZID" ), OrgDb = org.Hs.eg.db ) # shift gene symbol to entrez id
enrich <- enrichKEGG(gene =degID$ENTREZID ,
                     organism = 'hsa', 
                     pvalueCutoff=1, qvalueCutoff=1) # use OrgDb fitted for different species
GeneRatio <- as.numeric(lapply(strsplit(enrich$GeneRatio,split="/"),function(x) as.numeric(x[1])/as.numeric(x[2])))
BgRatio <- as.numeric(lapply(strsplit(enrich$BgRatio,split="/"),function(x) as.numeric(x[1])/as.numeric(x[2])  ))
enrich_factor <- GeneRatio/BgRatio
out <- data.frame(enrich$ID,enrich$Description,enrich$GeneRatio,enrich$BgRatio,round(enrich_factor,2),enrich$pvalue, enrich$p.adjust,enrich$qvalue,enrich$geneID)
colnames(out) <- c("ID","Description","GeneRatio","BgRatio","enrich_factor","pvalue", "p.adjust","qvalue","geneID")
setwd("/home/shengli/projects/AEG_proteomics/results/proteome")
write.table(out, "KEGG_Prot_diff_up.txt", row.names = F, sep="\t", quote = F)
### enrichment analysis for down-regulated proteins
rm(list = ls())
options(stringsAsFactors = FALSE) # prohibit shift form character to factor
library(tidyverse)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
dep_df <- read.table('Normal_Tumor_DEP_sig_genes.txt',header=T,row.names=1,sep='\t')
down_df <- dep_df[which(dep_df[,'logFC']<0),]
degID <- bitr(down_df$Gene, fromType = "SYMBOL", toType = c( "ENTREZID" ), OrgDb = org.Hs.eg.db ) # shift gene symbol to entrez id
enrich <- enrichKEGG(gene =degID$ENTREZID ,
                     organism = 'hsa', 
                     pvalueCutoff=1, qvalueCutoff=1) # use OrgDb fitted for different species
GeneRatio <- as.numeric(lapply(strsplit(enrich$GeneRatio,split="/"),function(x) as.numeric(x[1])/as.numeric(x[2])))
BgRatio <- as.numeric(lapply(strsplit(enrich$BgRatio,split="/"),function(x) as.numeric(x[1])/as.numeric(x[2])  ))
enrich_factor <- GeneRatio/BgRatio
out <- data.frame(enrich$ID,enrich$Description,enrich$GeneRatio,enrich$BgRatio,round(enrich_factor,2),enrich$pvalue, enrich$p.adjust,enrich$qvalue,enrich$geneID)
colnames(out) <- c("ID","Description","GeneRatio","BgRatio","enrich_factor","pvalue", "p.adjust","qvalue","geneID")
setwd("/home/shengli/projects/AEG_proteomics/results/proteome")
write.table(out, "KEGG_Prot_diff_down.txt", row.names = F, sep="\t", quote = F)        
## calculate hallmark scores at protein-level
rm(list=ls())
library(GSVA)
dir_res <- '/home/shengli/projects/AEG_proteomics/results/proteome'
hallmarks <- read.table('/home/public/public_data/MSigDB/processed/h.all.v6.2.symbols.txt',header=T,sep='\t',row.names=1)
dir_prot <- '/home/shengli/projects/AEG_proteomics/results/proteome'
file_aeg <- paste(dir_res,'/','Proteomics_iBAQ103_log2quantile_normlization_impute_perc25_gene.txt',sep='')
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
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
write.table(es_aeg,file='AEG_protein_hallmark_scores.txt',quote=F,sep='\t')

# difference of hallmark scores between tumor and normal samples
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
hk_scores <- read.table('AEG_protein_hallmark_scores.txt',header=T,row.names=1,sep='\t')
results_pvalues <- c()
results_foldchanges <- c()
for (n in 1:nrow(hk_scores)) {
  hk_tumor <- as.numeric(hk_scores[n,1:103])
  hk_normal <- as.numeric(hk_scores[n,104:206])
  mid_tumor <- median(hk_tumor)
  mid_normal <- median(hk_normal)
  fc <- round(mid_tumor/mid_normal,3)
  results_foldchanges <- c(results_foldchanges,fc)
  dif <- wilcox.test(hk_tumor,hk_normal)
  results_pvalues <- c(results_pvalues,dif$p.value)
}
results_fdr <- p.adjust(results_pvalues,method='fdr')
results <- cbind(results_foldchanges,results_pvalues,results_fdr)
rownames(results) <- rownames(hk_scores)
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
write.table(results,file='AEG_hallmark_diff.txt',sep='\t',quote=F,col.names=F)

### survival analysis for DEP

rm(list=ls())
library(survival)
aeg_clinical <- read.table('/home/shengli/projects/AEG_proteomics/data/clinical/survival_data_v3.txt',header=T,sep='\t')
dep_prot <- read.table('/home/shengli/projects/AEG_proteomics/results/proteome/Normal_Tumor_DEP_sig_genes.txt',header=T,sep='\t')
dep_prots <- dep_prot[,'Protein']
prot_profile <- read.table('/home/shengli/projects/AEG_proteomics/results/proteome/Proteomics_iBAQ103_log2quantile_normlization_impute_perc25.txt',header=T,row.names=1,sep='\t')
prot_profile <- prot_profile[,1:103]
rownames(aeg_clinical) <- aeg_clinical[,'Patient']
dep_profile <- prot_profile[dep_prots,]
aeg_clinical <- aeg_clinical[colnames(dep_profile),]
Results <- matrix(NA,nrow(dep_profile),11,byrow=TRUE)
rownames(Results) <- rownames(dep_profile)
colnames(Results) <- c('coef','Exp(coef)','Coxp','KMp','N','low','high','log2expr_FC','CI_95%_for_HR','Time_dif','log2Time_FC')

for (r in 1:nrow(Results)) {
  prot_exp <- as.numeric(dep_profile[r,])
  time_month <- aeg_clinical[,'Time']
  cen_status <- aeg_clinical[,'Status']
  Results[r,'N'] <- length(prot_exp)
  test_data1 <- list(time = time_month,
                     status = cen_status,
                     prot_expr = prot_exp)
  tempresult <- try(model0 <- coxph(Surv(time,status) ~ prot_expr, data=test_data1, na.action=na.exclude),silent=TRUE)
  if (!is(tempresult,'try-error')) {
    Results[r,c('coef','Exp(coef)','Coxp')] <- signif(as.numeric(summary(model0)$coefficients[1,c('coef','exp(coef)','Pr(>|z|)')]),digit=3)
    HR_detail = summary(model0)
    Results[r,c('CI_95%_for_HR')] = paste(signif(HR_detail$conf.int[,'lower .95'],4),'-',signif(HR_detail$conf.int[,'upper .95'],4),sep=' ')
    cutgroup <- ifelse(prot_exp <= median(prot_exp),'Low','High')
    if (length(unique(cutgroup)) > 1) {
      test_data2 <- list(time = time_month,
                         status = cen_status,
                         group = as.factor(cutgroup))
      model1 <- survdiff(Surv(time,status) ~ group, data=test_data2,na.action=na.exclude)
      Results[r,'KMp'] <- 1 - pchisq(model1$chisq,df=length(levels(factor(cutgroup))) - 1)
      Results[r,'low'] <- median(prot_exp[cutgroup == 'Low'])
      Results[r,'high'] <- median(prot_exp[cutgroup == 'High'])
      Results[r,'log2expr_FC'] <- log2((median(prot_exp[cutgroup == 'High'])+1)/(median(prot_exp[cutgroup == 'Low']) + 1))
      Results[r,'Time_dif'] <- mean(time_month[cutgroup == 'High']) - mean(time_month[cutgroup == 'Low'])
      Results[r,'log2Time_FC'] <- log2((median(time_month[cutgroup == 'High'])+1)/(median(time_month[cutgroup == 'Low']) + 1))
    }
  } else {
    cat("Coxph test doesn't pass!\n")
  }
}

tmp <- as.numeric(as.vector(Results[,'Coxp']))
CoxFDRAdjustPvalue <- p.adjust(tmp,method='fdr')
CoxOtherAdjustPvalue <- p.adjust(tmp,method='bonferroni')
tmp <- as.numeric(as.vector(Results[,'KMp']))
KMpFDRAdjustPvalue <- p.adjust(tmp,method='fdr')
KMpOtherAdjustPvalue <- p.adjust(tmp,method='bonferroni')
Results <- cbind(Marker=rownames(Results),Results,FDR=CoxFDRAdjustPvalue,Bonferroni=CoxOtherAdjustPvalue,
                KMFDR=KMpFDRAdjustPvalue, KMBonferroni=KMpOtherAdjustPvalue)

write.table(Results,file='/home/shengli/projects/AEG_proteomics/results/proteome/survival_analysis/DEP_survival.txt',sep='\t',quote=FALSE,row.names=FALSE)


