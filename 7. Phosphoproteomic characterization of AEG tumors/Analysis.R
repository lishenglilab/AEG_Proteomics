### normalization
library(limma)
setwd('/home/shengli/projects/AEG_proteomics/data/phosphoproteome')
prot_ibaq <- read.table('Phosphoproteomics_iBAQ103.txt',header=T,row.names=1)
prot_ibaq <- as.matrix(prot_ibaq)
quant_norm_ibaq <- normalizeQuantiles(prot_ibaq)
write.table(quant_norm_ibaq,file='Phosphoproteomics_iBAQ103_quantile_normlization.txt',quote=F,sep='\t')
log2_norm_ibaq <- log2(quant_norm_ibaq)
write.table(log2_norm_ibaq,file='Phosphoproteomics_iBAQ103_log2quantile_normlization.txt',quote=F,sep='\t')

## differential phosphorylate sites
library(limma)
setwd('/home/shengli/projects/AEG_proteomics/data/phosphoproteome')
phosp_mx <- read.table('Phosphoproteomics_iBAQ103_log2quantile_normlization_impute_perc25.txt',header=T,row.names=1,sep='\t')
design_nt <- cbind(Intercept=1,
                   Group=c(rep(1,103),rep(0,103)))
fit <- lmFit(phosp_mx,design_nt)
fit <- eBayes(fit,trend=TRUE)
dep_nt <- topTable(fit,coef=2,number=Inf)
setwd('/home/shengli/projects/AEG_proteomics/results/phosphoproteome')
write.table(dep_nt,file="Normal_Tumor_DEPhop_all.txt",quote=F,sep='\t',append=F)

## differential phosphoprotomics expression in each subtype between tumor and NAT samples
library(limma)
setwd('/home/shengli/projects/AEG_proteomics/data/phosphoproteome')
phosp_mx <- read.table('Phosphoproteomics_iBAQ103_log2quantile_normlization_impute_perc25.txt',header=T,row.names=1,sep='\t')
samples_phosp <- colnames(phosp_mx)
#S-I subtype
samples_t1 <- c('C0551594','C0559459','C0560016','C0563560','C0564447','C0569091','C0575352','C0579894','C0582702','C0584162','C0586163','C0583851','C0590054','C0595498',
                'C0596146','C1005614','C1006005','C1009828','C0532443','C0505919','C0562401','C0587488','C0599521','C0584104','C0498355','C0565023','C246567',
                'C0576092','C0488963','C0533179','C0504861','C0568971','C0528365','C0533153','C1007477','C0572176','C0594667','C0522528','C0561723','C0584243')
samples_t1_p <- samples_t1[which(samples_t1 %in% samples_phosp)]
samples_n1 <- c('N0551594','N0559459','N0560016','N0563560','N0564447','N0569091','N0575352','N0579894','N0582702','N0584162','N0586163','N0583851','N0590054','N0595498',
                'N0596146','N1005614','N1006005','N1009828','N0532443','N0505919','N0562401','N0587488','N0599521','N0584104','N0498355','N0565023','N246567',
                'N0576092','N0488963','N0533179','N0504861','N0568971','N0528365','N0533153','N1007477','N0572176','N0594667','N0522528','N0561723','N0584243')
samples_n1_p <- samples_n1[which(samples_n1 %in% samples_phosp)]

mx_s1 <- phosp_mx[,c(samples_t1_p,samples_n1_p)]
design_s1 <- cbind(Intercept=1,
                   Group=c(rep(1,length(samples_t1_p)),rep(0,length(samples_t1_p))))
fit <- lmFit(mx_s1,design_s1)
fit <- eBayes(fit,trend = TRUE)
deg_s1 <- topTable(fit,coef=2,number=Inf)
setwd('/home/shengli/projects/AEG_proteomics/results/phosphoproteome')
write.table(deg_s1,file='Phosp_S1_diff_all.txt',quote=F,sep='\t',append=F)

# enrichment analysis of differential phosphorylation sites
rm(list = ls())
options(stringsAsFactors = FALSE) # prohibit shift form character to factor
library(tidyverse)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd('/home/shengli/projects/AEG_proteomics/results/phosphoproteome')
deg_df <- read.table('Phosp_S1_diff_sig_gene.txt',head=T,row.names=1,sep='\t')
degID <- bitr(deg_df$Gene, fromType = "SYMBOL", toType = c( "ENTREZID" ), OrgDb = org.Hs.eg.db ) # shift gene symbol to entrez id
enrich <- enrichKEGG(gene =degID$ENTREZID ,
                     organism = 'hsa', 
                     pvalueCutoff=1, qvalueCutoff=1) # use OrgDb fitted for different species
GeneRatio <- as.numeric(lapply(strsplit(enrich$GeneRatio,split="/"),function(x) as.numeric(x[1])/as.numeric(x[2])))
BgRatio <- as.numeric(lapply(strsplit(enrich$BgRatio,split="/"),function(x) as.numeric(x[1])/as.numeric(x[2])  ))
enrich_factor <- GeneRatio/BgRatio
out <- data.frame(enrich$ID,enrich$Description,enrich$GeneRatio,enrich$BgRatio,round(enrich_factor,2),enrich$pvalue, enrich$p.adjust,enrich$qvalue,enrich$geneID)
colnames(out) <- c("ID","Description","GeneRatio","BgRatio","enrich_factor","pvalue", "p.adjust","qvalue","geneID")
setwd("/home/shengli/projects/AEG_proteomics/results/phosphoproteome")
write.table(out, "KEGG_Phosp_S1_diff.txt", row.names = F, sep="\t", quote = F)

## KSEA analysis
library(KSEAapp)
setwd('/home/shengli/projects/AEG_proteomics/results/phosphoproteome')
px <- read.table('Phosp_S3_diff_sig_gene_KSEA.csv',header=T,sep=',')
KSEA.Complete(KSData,PX=px,NetworKIN = T,NetworKIN.cutoff = 1,m.cutoff=5,p.cutoff = 0.05)

## calculate correlations between kinases and substrates
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
prot_mx <- read.table('Proteomics_iBAQ_log2quantile_normlization_impute_tumor_perc25_gene.txt',header=T,row.names=1,sep='\t')
setwd('/home/shengli/projects/AEG_proteomics/data/phosphoproteome')
phosph_mx <- read.table('Phosphoproteomics_iBAQ103_log2quantile_normlization_impute_perc25.txt',header=T,row.names=1,sep='\t')

#S-I subtype
samples_t1 <- c('C0551594','C0559459','C0560016','C0563560','C0564447','C0569091','C0575352','C0579894','C0582702','C0584162','C0586163','C0583851','C0590054','C0595498',
                'C0596146','C1005614','C1006005','C1009828','C0532443','C0505919','C0562401','C0587488','C0599521','C0584104','C0498355','C0565023','C246567',
                'C0576092','C0488963','C0533179','C0504861','C0568971','C0528365','C0533153','C1007477','C0572176','C0594667','C0522528','C0561723','C0584243')
setwd('/home/shengli/projects/AEG_proteomics/results/phosphoproteome')
links_t1 <- read.table('Kinase_substrate_sig_S1.txt',header=T,sep='\t')
kinases_t1 <- unique(links_t1[,'Kinase'])
substrate_t1 <- unique(links_t1[,'Substrate_2'])
prot_mx_t1 <- prot_mx[kinases_t1,samples_t1]
phosph_mx_t1 <- phosph_mx[substrate_t1,samples_t1]

Results <- matrix(NA,nrow=nrow(links_t1),ncol=4,byrow=TRUE)
colnames(Results) <- c('Kinase','Substrate','Correlation','Pvalue')
for (r in 1:nrow(links_t1)) {
  cat(r,'\n')
  kin <- links_t1[r,'Kinase']
  subs <- links_t1[r,'Substrate_2']
  expr_prot <- as.numeric(prot_mx_t1[kin,])
  expr_phosph <- as.numeric(phosph_mx_t1[subs,])
  if (is.na(expr_prot) | is.na(expr_phosph)) {
    Results[r,] <- NA
    next
  }
  test_expr <- cor.test(expr_prot,expr_phosph,method='spearman')
  Results[r,'Kinase'] <- kin
  Results[r,'Substrate'] <- subs
  Results[r,'Correlation'] <- test_expr$estimate
  Results[r,'Pvalue'] <- test_expr$p.value
}

setwd('/home/shengli/projects/0.Collaborations/4.Qin.Jiangjiang/results/phosphoproteome')
write.table(Results,file='Kinase_substrate_cor_S1.txt',sep='\t',row.names=F,quote=F)


# S-II
samples_t2 <- c('C217360','C266198','C325352','C227213','C0502774','C284805','C0538376','C335689','C242946','C242956','C0505805','C194922','C197023','C195615',
                'C225468','C201779','C260052','C265603','C0503345','C0502439','C0540265','C234760','C219685')
setwd('/home/shengli/projects/AEG_proteomics/results/phosphoproteome')
links_t2 <- read.table('Kinase_substrate_sig_S2.txt',header=T,sep='\t')
kinases_t2 <- unique(links_t2[,'Kinase'])
substrate_t2 <- unique(links_t2[,'Substrate_2'])
prot_mx_t2 <- prot_mx[kinases_t2,samples_t2]
phosph_mx_t2 <- phosph_mx[substrate_t2,samples_t2]

Results <- matrix(NA,nrow=nrow(links_t2),ncol=4,byrow=TRUE)
colnames(Results) <- c('Kinase','Substrate','Correlation','Pvalue')
for (r in 1:nrow(links_t2)) {
  cat(r,'\n')
  kin <- links_t2[r,'Kinase']
  subs <- links_t2[r,'Substrate_2']
  expr_prot <- as.numeric(prot_mx_t2[kin,])
  expr_phosph <- as.numeric(phosph_mx_t2[subs,])
  if (is.na(expr_prot) | is.na(expr_phosph)) {
    Results[r,] <- NA
    next
  }
  test_expr <- cor.test(expr_prot,expr_phosph,method='spearman')
  Results[r,'Kinase'] <- kin
  Results[r,'Substrate'] <- subs
  Results[r,'Correlation'] <- test_expr$estimate
  Results[r,'Pvalue'] <- test_expr$p.value
}

setwd('/home/shengli/projects/AEG_proteomics/results/phosphoproteome')
write.table(Results,file='Kinase_substrate_cor_S2.txt',sep='\t',row.names=F,quote=F)

# S-III subtype
samples_t3 <- c('C332654','C386485','C266473','C338321','C338421','C335207','C378018','C325830','C196713','C341159','C151058','C203659','C203933','C269944',
                'C210038','C260824','C218213','C391965','C363245','C196436','C321599','C241209','C215502','C248148','C322712','C248917','C0600290','C291416',
                'C1004644','C305692','C258083','C259302','C1003009','C288415','C314628','C260319','C199181','C256818','C260513','C260742')
setwd('/home/shengli/projects/AEG_proteomics/results/phosphoproteome')
links_t3 <- read.table('Kinase_substrate_sig_S3.txt',header=T,sep='\t')
kinases_t3 <- unique(links_t3[,'Kinase'])
substrate_t3 <- unique(links_t3[,'Substrate_2'])
prot_mx_t3 <- prot_mx[kinases_t3,samples_t3]
phosph_mx_t3 <- phosph_mx[substrate_t3,samples_t3]

Results <- matrix(NA,nrow=nrow(links_t3),ncol=4,byrow=TRUE)
colnames(Results) <- c('Kinase','Substrate','Correlation','Pvalue')
for (r in 1:nrow(links_t3)) {
  cat(r,'\n')
  kin <- links_t3[r,'Kinase']
  subs <- links_t3[r,'Substrate_2']
  expr_prot <- as.numeric(prot_mx_t3[kin,])
  expr_phosph <- as.numeric(phosph_mx_t3[subs,])
  if (is.na(expr_prot) | is.na(expr_phosph)) {
    Results[r,] <- NA
    next
  }
  test_expr <- cor.test(expr_prot,expr_phosph,method='spearman')
  Results[r,'Kinase'] <- kin
  Results[r,'Substrate'] <- subs
  Results[r,'Correlation'] <- test_expr$estimate
  Results[r,'Pvalue'] <- test_expr$p.value
}

setwd('/home/shengli/projects/AEG_proteomics/results/phosphoproteome')
write.table(Results,file='Kinase_substrate_cor_S3.txt',sep='\t',row.names=F,quote=F)
                

