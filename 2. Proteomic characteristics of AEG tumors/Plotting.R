### Volcano plot showing the difference in proteins between AEG tumor and paired NAT samples. Fig. 2a

## valcono plot of differentially expressed proteins, Fig. 2a
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

# heatmap for hallmark scores
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
hk_scores <- read.table('AEG_protein_hallmark_scores.txt',header=T,row.names=1,sep='\t')
pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics/Prot_hk_scores_heatmap.pdf',width=5,height=8)
pheatmap(hk_scores,scale='row',show_rownames=T,show_colnames=F,cluster_col=F,cluster_row=T,fontsize_col=5,fontsize_row=6,border_color=NA,color=colorRampPalette(c('blue','white','red'))(100))
dev.off()

# boxplot for individual hallmarks, Fig. 2d
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

## plot survival curves for individual hallmark scores, 
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

### heatmap for top DEPs in PPI network, Fig. 2f
prots <- c('P42771','P10145','P62308','P11215','P38919','P17844','Q92841','P07910','P12004','P09651','P46531','P05067','P08238',
           'P42224','P20226','Q07817','O14980','P49915','P35354','P63165','P04637','P07900','P02751','Q13547','P01730','Q09472',
           'O43143','P06493','P22087','P40763','P14635','Q07955','Q15029','P11388','Q15910','Q96EB6','Q9BZE4','P23921','Q9Y265',
           'P24941','P19338','P14780','P04626')
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
prot_mx <- read.table('Proteomics_iBAQ103_log2quantile_normlization_impute_perc25.txt',header=T,row.names=1,sep='\t')
prot_mx <- prot_mx[prots,]
prot_mx <- as.matrix(prot_mx)
pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics/DEP_PPI_heatmap.pdf',width=5,height=8)
pheatmap(prot_mx,scale='row',show_rownames=T,show_colnames=F,cluster_col=F,cluster_row=T,fontsize_col=5,fontsize_row=10,border_color=NA,color=colorRampPalette(c('blue','white','red'))(100))
dev.off()

### heatmap for differential druggable genes, Supplementary Fig. 5a
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
dep_drg <- read.table('DEP_sig_druggable_genes.txt',header=T)
dep_drg2 <- dep_drg[which(abs(dep_drg[,'logFC']) > 1),]
prots <- c(rownames(dep_sig_mx)[1:49],'P35869')
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
prot_mx <- read.table('Proteomics_iBAQ103_log2quantile_normlization_impute_perc25.txt',header=T,row.names=1,sep='\t')

dep_sig_mx <- prot_mx[prots,]

dep_sig_mx <- as.matrix(dep_sig_mx)

pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics/DEP_druggable_heatmap_v2.pdf',width=5,height=8)
pheatmap(dep_sig_mx,scale='row',show_rownames=T,show_colnames=F,cluster_col=F,cluster_row=T,fontsize_col=5,fontsize_row=6,border_color=NA,color=colorRampPalette(c('blue','white','red'))(100))
dev.off()
### boxplot for individual proteins, Supplementary Fig. 5b
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
prot_mx <- read.table('Proteomics_iBAQ103_log2quantile_normlization_impute_perc25.txt',header=T,row.names=1,sep='\t')
prot <- 'P35869'
expr_normal <- as.numeric(prot_mx[prot,104:206])
expr_tumor <- as.numeric(prot_mx[prot,1:103])
expr <- data.frame(Expr=c(expr_normal,expr_tumor),
                     Sample=c(colnames(prot_mx)[104:206],colnames(prot_mx)[1:103]),
                     Group=c(rep('Normal',103),rep('Tumor',103)))
pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics/Proteomics_P35869_nt.pdf',width=5,height = 4)
ggplot(expr,aes(x=Group,y=Expr)) +
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



