### Volcano plot showing the difference in proteins between AEG tumor and paired NAT samples. Fig. 2a

## valcono plot of differentially expressed proteins
rm(list=ls())
library(ggplot2)
library(ggrepel)
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
diff <- read.table('Normal_Tumor_DEP_all.txt',header=T,row.names=1,sep='\t')
diff$lgfdr <- -log10(diff[,'adj.P.Val'])
diff$prot <- row.names(diff)
diff_none <- diff[which(diff[,'adj.P.Val'] >= 0.01 | abs(diff[,'logFC']) <= 1),]
diff_up <- diff[which(diff[,'adj.P.Val'] < 0.01 & diff[,'logFC'] > 1),]
diff_down <- diff[which(diff[,'adj.P.Val'] < 0.01 & diff[,'logFC'] < -1),]
up_prot <- as.character(row.names(diff_up))
down_prot <- as.character(row.names(diff_down))

pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics/diff_prot_volcano_R1.pdf',width=6.87,height=5.92)
ggplot(diff,aes(x=logFC,y=lgfdr,label=prot)) +
  geom_point(data=diff_none,aes(x=logFC,y=lgfdr),color='lightgray') +
  geom_point(data=diff_up,aes(x=logFC,y=lgfdr),color='red') +
  geom_point(data=diff_down,aes(x=logFC,y=lgfdr),color='blue') +
  geom_label_repel(data = subset(diff, lgfdr > 30 | abs(logFC) > 8),
                   size = 3,
                   box.padding = 0.5,
                   point.padding = 0.5,
                   segment.size = 0.2,
                   segment.color = "grey50") +
  theme(panel.background=element_rect(colour="black",fill="white"))
dev.off()


### Functional enrichment results of up-regulated and down-regulated proteins, respectively. Fig. 2b
## enrichment analysis for up-regulated proteins
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
## barplot 
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
res_enrich <- read.table('KEGG_Prot_diff_up.txt',header=T,sep='\t')
enrich_4bar <- data.frame(bp=c('Spliceosome','DNA replication','Ribosome biogenesis in eukaryotes','Cell cycle','Nucleotide excision repair',
                               'Basal transcription factors','Nucleocytoplasmic transport','Mismatch repair','Phagosome','Human T-cell leukemia virus 1 infection'),
                          fdr=c(4.15809612261286e-21,1.41577604969151e-15,2.88839231289716e-08,5.67494339107914e-08,6.78767919267188e-08,
                                1.34148845574837e-07,1.43595014672821e-07,1.86237820286443e-07,1.19333339307256e-06,1.61648208424137e-05))
enrich_4bar$lgp <- -log10(enrich_4bar$fdr)
rownames(enrich_4bar) <- enrich_4bar[,'bp']
pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics/diff_up_prot_KEGG.pdf',width=7,height = 5)
barplot(enrich_4bar[,'lgp'],horiz=T,names.arg=rownames(enrich_4bar),las=2,xlim=c(0,20))
dev.off()

## down-regulated proteins
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

## barplot
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
res_enrich <- read.table('KEGG_Prot_diff_down.txt',header=T,sep='\t')
enrich_4bar <- data.frame(bp=c('Oxidative phosphorylation','Chemical carcinogenesis - reactive oxygen species','Valine, leucine and isoleucine degradation',
                               'Carbon metabolism','Propanoate metabolism','Retrograde endocannabinoid signaling','Citrate cycle (TCA cycle)','Butanoate metabolism',
                               'Glyoxylate and dicarboxylate metabolism','Pyruvate metabolism'),
                          fdr=c(1.25249915786254e-43,1.02713530477128e-32,1.41157794076512e-19,1.21780142137805e-17,6.03883251719555e-14,6.03883251719555e-14,
                                2.47024346542946e-12,1.36937413955737e-09,6.60964456641259e-08,1.10916540941525e-07))
enrich_4bar$lgp <- -log10(enrich_4bar$fdr)
rownames(enrich_4bar) <- enrich_4bar[,'bp']
pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics/diff_down_prot_KEGG.pdf',width=7,height = 5)
barplot(enrich_4bar[,'lgp'],horiz=T,names.arg=rownames(enrich_4bar),las=2,xlim=c(0,45))
dev.off()
                             
### heatmap of hallmark scores  Fig. 2c
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

# heatmap for hallmark scores
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
hk_scores <- read.table('AEG_protein_hallmark_scores.txt',header=T,row.names=1,sep='\t')
pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics/Prot_hk_scores_heatmap.pdf',width=5,height=8)
pheatmap(hk_scores,scale='row',show_rownames=T,show_colnames=F,cluster_col=F,cluster_row=T,fontsize_col=5,fontsize_row=6,border_color=NA,color=colorRampPalette(c('blue','white','red'))(100))
dev.off()

# barplot for individual hallmarks
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
hk_scores <- read.table('AEG_protein_hallmark_scores.txt',header=T,row.names=1,sep='\t')
hk <- 'HALLMARK_APICAL_JUNCTION'
scores_normal <- as.numeric(hk_scores[hk,104:206])
scores_tumor <- as.numeric(hk_scores[hk,1:103])
scores <- data.frame(Score=c(scores_normal,scores_tumor),
                     Sample=c(colnames(hk_scores)[104:206],colnames(hk_scores)[1:103]),
                     Group=c(rep('Normal',103),rep('Tumor',103)))
pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics/hk_apical_junction_nt.pdf',width=5,height = 4)
ggplot(scores,aes(x=Group,y=Score)) +
  geom_jitter(width=0.15,size=1.5) +
  geom_boxplot(color="black",fill=NA,outlier.shape=NA,width=0.3) +
  scale_x_discrete(limit=c("Normal","Tumor"),expand=c(0.1,0)) +
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

## plot survival curves for individual hallmark scores
library(survival)
library(ggplot2)
library(survminer)
setwd('/home/shengli/projects/AEG_proteomics/data/clinical')
surv_data <- read.table('survival_data_v3.txt',header=T,sep='\t')
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
hk_scoress <- read.table('AEG_protein_hallmark_scores.txt',header=T,row.names=1,sep='\t')
hk_sc_tumor <- hk_scoress['HALLMARK_KRAS_SIGNALING_UP',1:103]
mid_sc <- median(as.numeric(hk_sc_tumor))
pt_high <- colnames(hk_sc_tumor)[which(hk_sc_tumor[1,] > mid_sc)]
pt_low <- colnames(hk_sc_tumor)[which(hk_sc_tumor[1,] <= mid_sc)]
surv_data[which(surv_data[,'Patient'] %in% pt_high),'Group'] <- 'High'
surv_data[which(surv_data[,'Patient'] %in% pt_low),'Group'] <- 'Low'
fit <- survfit(Surv(Time,Status)~Group, data=surv_data)
 
pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics/hk_kras_signaling_up_survival.pdf',width=5.33,height=4.82)
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

