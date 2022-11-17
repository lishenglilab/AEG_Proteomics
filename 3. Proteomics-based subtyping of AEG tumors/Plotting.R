# heatmap for subtypes DEPs, Fig. 3a
rm(list=ls())
library(pheatmap)
setwd('/home/shengli/projects/AEG_proteomics/results')
prn_profile <- read.table('Proteomics_perc25_impute_cvs25_v2.txt',header=T,row.names=1,sep='\t')
setwd('/home/shengli/projects/AEG_proteomics/data/clinical')
subtypes <- read.table('survival_data_v3.txt',header=T,row.names=1,sep='\t')
samples_c1 <- as.character(rownames(subtypes[which(subtypes[,'Group'] == 'Group1'),]))
samples_c2 <- as.character(rownames(subtypes[which(subtypes[,'Group'] == 'Group2'),]))
samples_c3 <- as.character(rownames(subtypes[which(subtypes[,'Group'] == 'Group3'),]))
samples_c1c2c3 <- c(samples_c1,samples_c2,samples_c3)
setwd('/home/shengli/projects/AEG_proteomics/results')
dep_c1c2 <- read.table('G1_G2_DEP_all.txt',header=T,row.names=1,sep='\t')
sigp_c1c2 <- as.character(rownames(dep_c1c2[which(abs(dep_c1c2[,'logFC'])>1 & dep_c1c2[,'adj.P.Val'] < 0.05),]))
dep_c2c3 <- read.table('G2_G3_DEP_all.txt',header=T,row.names=1,sep='\t')
sigp_c2c3 <- as.character(rownames(dep_c2c3[which(abs(dep_c2c3[,'logFC'])>1 & dep_c2c3[,'adj.P.Val'] < 0.05),]))
dep_c1c3 <- read.table('G1_G3_DEP_all.txt',header=T,row.names=1,sep='\t')
sigp_c1c3 <- as.character(rownames(dep_c1c3[which(abs(dep_c1c3[,'logFC'])>1 & dep_c1c3[,'adj.P.Val'] < 0.05),]))
sigps <- c(sigp_c1c2,sigp_c2c3,sigp_c1c3)
sigps <- unique(sigps)
prn_order <- prn_profile[sigps,samples_c1c2c3]
pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics_subtypes_heatmap.pdf',height=5.5, width=6.5)
pheatmap(prn_order,scale='row',show_rownames=F,show_colnames=T,cluster_col=F,cluster_row=T,color=colorRampPalette(c('blue','white','red'))(100))
dev.off()

## plot individual survival curves of proteins, Supplementary Fig. 10
library(survival)
library(ggplot2)
library(survminer)
setwd('/home/shengli/projects/AEG_protoemics/results/proteome')
data_sigp <- read.table('sigp_for_coxp_v2.txt',header=T,sep='\t')

protein <- 'Q96AT1'
abundances <- as.numeric(data_sigp[,protein])
median_abund <- median(abundances)
data_sigp$Group <- rep('High',nrow(data_sigp))
data_sigp[which(data_sigp[,protein] <= median_abund),'Group'] <- 'Low'
fit <- survfit(Surv(Time,Status)~Group, data=data_sigp)
cox_model <- coxph(Surv(Time,Status)~Group+Age+Sex+Smoking+Alcohol+Siewert_type+Clinical_stage, data=data_sigp)
summary(cox_model)
# summary(cox_model) HR=0.48 (95% CI, 0.26-0.88)
pdf('/home/shengli/projects/AEG_proteomics/figures/Survival_curve_Q96AT1.pdf',height=5,width=5)
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

#plot all difference of gene-level hallmarks, Fig. 3d
setwd('/home/shengli/projects/AEG_proteomics/results/RNAseq')
hk_dif <- read.table('AEG_hallmarks_diff.txt',header=T,sep='\t')
hallmarks <- c('Apical junction','Apical surface','Peroxisome',
               'Adipogenesis','Angiogenesis','EMT','Myogenesis','Spermatogenesis','Pancreas beta cells',
               'DNA repair','UV response UP','UV response DN',
               'Allograft rejection','Coagulation','Complement','Interferon alpha response','Interferon gamma response','IL6 JAK STAT3 signaling','Inflammatory response',
               'Bile acid metabolism','Cholesterol homeostasis','Fatty acid metabolism','Glycolysis','HEME metabolism','Oxidative phosphorylation','Xenobiotic metabolism',
               'Apoptosis','Hypoxia','Protein secretion','Unfolded protein response','Reactive oxigen species pathway',
               'E2F targets','G2M checkpoint','MYC targets v1','MYC targets v2','P53 pathway','Mitotic spindle',
               'Androgen response','Estrogen response early','Estrogen response late','IL2 STAT5 signaling','KRAS signaling UP','KRAS signaling DN','mTORC1 signaling','Notch signaling',
               'PI3K AKT mTOR signaling','Hedgehog signaling','TGF beta signaling','TNFA signaling via NFKB','WNT beta catenin signaling')
hk_dif$lgfdr <- -log10(hk_dif[,'FDR'])
hk_dif$logfc <- log2(hk_dif[,'FoldChange'])
pdf('/home/shengli/projects/AEG_proteomics/figures/AEG_hk_scores_diff.pdf',width=9.5,height=3.5)
ggplot(hk_dif,aes(x=Hallmark,y=Subtype))+
  geom_point(aes(size=lgfdr,col=logfc))+
  scale_color_gradient2(low="blue",mid="white",high="red",midpoint=0,na.value="white",name="FoldChange")+
  scale_size_continuous(limit=c(0,15),range=c(0.5,4),breaks=c(-log10(0.05),5,10,15),labels=c("0.05","1e-5","1e-10","<1e-15"))+
  scale_y_discrete(limit=c('S1','S2','S3'),expand=c(0.05,0.05))+
  scale_x_discrete(limit=hallmarks,expand=c(0.05,0.05)) +
  theme(panel.background=element_rect(colour="black",fill="white"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=11,colour="black"),
        axis.text.x=element_text(size=10,colour="black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.key=element_rect(fill="white",colour="black"))
dev.off()

# plot single hallmarks, Fig. 3e-f
setwd('/home/shengli/projects/AEG_proteomics/results/RNAseq')
hallmark_scores <- read.table('AEG_hallmark_scores.txt',header=T,row.names=1,sep='\t')
scores_s1 <- as.numeric(hallmark_scores['HALLMARK_PANCREAS_BETA_CELLS',samples_t1_p])
scores_s2 <- as.numeric(hallmark_scores['HALLMARK_PANCREAS_BETA_CELLS',samples_t2_p])
scores_s3 <- as.numeric(hallmark_scores['HALLMARK_PANCREAS_BETA_CELLS',samples_t3_p])

scores <- data.frame(Score=c(scores_s1,scores_s2,scores_s3),
                     Sample=c(samples_t1_p,samples_t2_p,samples_t3_p),
                     Group=c(rep('S1',38),rep('S2',16),rep('S3',30)))
pdf('/home/shengli/projects/AEG_proteomics/figures/RNAseq/hk_pancreas_beta_cells_comp.pdf',width=5,height = 4)
ggplot(scores,aes(x=Group,y=Score)) +
  geom_jitter(width=0.15,size=1.5) +
  geom_boxplot(color="black",fill=NA,outlier.shape=NA,width=0.3) +
  scale_x_discrete(limit=c("S1","S2","S3"),expand=c(0.1,0)) +
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

## plot multi-variate Cox regression analysis of clinical features, Supplementary Fig. 7
setwd('/home/shengli/projects/AEG_proteomcs/data/clinical')
data_clinical <- read.table('survival_data_v4.txt',header=T,sep='\t')
model_cox <- coxph(Surv(Time,Status) ~ Group + Age + Sex + Smoking + Alcohol + Siewert_type + Clinical_stage, data = data_clinical)
pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics/clinical_multicox.pdf',height=3.5,width=6.5)
ggforest(model_cox)
dev.off()

### heatmap for signature proteins, Fig. 3g
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
sig_prf <- read.table('Proteomics_signature_proteins.txt',header=T,row.names=1,sep='\t')
res <- read.table('signature_proteins_cox.txt',header=T,row.names=1,sep='\t')
sigp_cox <- res[which(res[,'p.value'] < 0.05),]

prt_order <- c('P35670','Q96GY3','Q92994','P02775','Q96AT1','Q9H4M3','P24592','Q6H8Q1','Q13563','Q6FHJ7','P29279','P04234')
pdf('/home/shengli/projects/AEG_proteomics/figures/sigp_coxsig_heatmap.pdf',height=3.5,width=7)
pheatmap(sig_prf[prt_order,],scale='row',show_rownames=T,show_colnames=T,cluster_col=F,cluster_row=F,fontsize_col=5,fontsize_row=5,color=colorRampPalette(c('blue','white','red'))(100))
dev.off()
### cox forest plot
headers <- c('Patient','Status','Time',prt_order)
data_sigp <- data_mx[,headers]
model_cox <- coxph(Surv(Time,Status) ~ P35670 + Q96GY3 + Q92994 + P02775 + Q96AT1 + Q9H4M3 + P24592 + Q6H8Q1 + Q13563 + Q6FHJ7 + P29279 + P04234, data = data_sigp)
pdf('/home/shengli/projects/AEG_proteomics/figures/sigp_multicox.pdf',height=3.5,width=6.5)
ggforest(model_cox)
dev.off()

## waterfall plot in each subtype, Supplementary Fig. 8
library(maftools)
setwd('/home/shengli/projects/AEG_proteomics/results/WES')
aeg_s3 <- read.maf(maf='AEG_mutation_neat_filter_non151058_group3.maf')
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/mutation_waterfall_S3_top30.pdf',height=6,width=10)
oncoplot(maf=aeg_s3,draw_titv=T,top=30)
dev.off()

# pie chart of mutation-protein affects, Supplementary Fig.9a
s1_effects <- c(24453,40731)
names(s1_effects) <- c('up','down')
pdf('/home/shengli/projects/AEG_proteomics/figures/2022.08/AEG_mut_prot_effect_stat_s1.pdf',height = 5, width =5)
pie(s1_effects)
dev.off()

s2_effects <- c(1328,2572)
names(s2_effects) <- c('up','down')
pdf('/home/shengli/projects/AEG_proteomics/figures/2022.08/AEG_mut_prot_effect_stat_s2.pdf',height=5,width=5)
pie(s2_effects)
dev.off()


s3_effects <- c(432,714)
names(s3_effects) <- c('up','down')
pdf('/home/shengli/projects/AEG_proteomics/figures/2022.08/AEG_mut_prot_effect_stat_s3.pdf',height = 5, width = 5)
pie(s3_effects)
dev.off()

