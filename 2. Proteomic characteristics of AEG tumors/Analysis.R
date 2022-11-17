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

