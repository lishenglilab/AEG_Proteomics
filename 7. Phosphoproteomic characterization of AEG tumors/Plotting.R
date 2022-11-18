### phosphorylation abundance distribution, Supplementary Fig. 1b
rm(list=ls())
library(ggplot2)
setwd('/home/shengli/projects/AEG_proteomics/data/phosphoproteome')
prot_ibaq <- read.table('Phosphoproteomics_iBAQ103_log2quantile_normlization_df.txt',header=T)
sample_cancer <- as.character(unique(prot_ibaq[which(prot_ibaq[,'Group']=='Cancer'),'Sample']))
sample_normal <- as.character(unique(prot_ibaq[which(prot_ibaq[,'Group']=='Normal'),'Sample']))
samples <- c(sample_cancer,sample_normal)
pdf('/home/shengli/projects/AEG_proteomics/figures/Phosphoproteomics_iBAQ_norm_distribution.pdf',width=20,height=5)
ggplot(prot_ibaq,aes(x=Sample,y=Abundance,fill=Group)) +
  geom_boxplot()
dev.off()

## phosphosrylation site numbers, Fig. 1d
setwd('/home/shengli/projects/AEG_proteomics/results/phosphoproteome')
phosp_num <- read.table('Phosphoproteomics_sample_numbers.txt',header=T,sep='\t')
num_normal <- phosp_num[which(phosp_num[,'Type']=='normal'),]
num_tumor <- phosp_num[which(phosp_num[,'Type']=='tumor'),]
rownames(num_normal) <- num_normal[,'Sample']
rownames(num_tumor) <- num_tumor[,'Sample']
num_normal_ordered <- num_normal[order(num_normal[,'Phosphoproteomics_number'],decreasing = F),]
num_tumor_ordered <- num_tumor[order(num_tumor[,'Phosphoproteomics_number'],decreasing = F),]
pdf('/home/shengli/projects/AEG_proteomics/figures/Phosphoproteomics/Phospho_numbers_tumor.pdf',height=5,width=5)
plot(num_tumor_ordered[,'Phosphoproteomics_number'],ylim=c(5500,10000),las=1)
dev.off()

num_normal <- phosp_num[which(phosp_num[,'Type']=='normal'),]
rownames(num_normal) <- num_normal[,'Sample']
num_normal_ordered <- num_normal[rownames(num_tumor_ordered),]
pdf('/home/shengli/projects/AEG_proteomics/figures/Phosphoproteomics/Protein_numbers_normal.pdf',height=5,width=5)
plot(num_normal_ordered[,'Protein_number'],ylim=c(8000,9500),las=1)
dev.off()

library(ggplot2)
nums <- c(as.numeric(num_tumor_ordered[,'Phosphoproteomics_number']),as.numeric(num_normal_ordered[,'Phosphoproteomics_number']))
nums_df <- data.frame(Numbers=nums,
                      sample_group=c(rep('tumor',103),rep('normal',103)))
pdf('/home/shengli/projects/AEG_proteomics/figures/Phosphoproteomics/phosp_numbers_comp.pdf',height=6,width=5)
ggplot(nums_df,aes(x=sample_group,y=Numbers)) +
  geom_jitter(width=0.15,size=0.7) +
  geom_boxplot(color="black",fill=NA,outlier.shape=NA,width=0.3) +
  scale_x_discrete(limit=c('normal','tumor'),expand=c(0.1,0)) +
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
## valcano plot of differential phosphorylation sites, Fig. 7a
rm(list=ls())
library(ggplot2)
library(ggrepel)
setwd('/home/shengli/projects/AEG_proteomics/results/phosphoproteome')
diff <- read.table('Normal_Tumor_DEPhop_all.txt',header=T,sep='\t')
diff$lgfdr <- -log10(diff[,'adj.P.Val'])
diff$phosp <- row.names(diff)
diff_none <- diff[which(diff[,'adj.P.Val'] >= 0.01 | abs(diff[,'logFC']) <= 1),]
diff_up <- diff[which(diff[,'adj.P.Val'] < 0.01 & diff[,'logFC'] > 1),]
diff_down <- diff[which(diff[,'adj.P.Val'] < 0.01 & diff[,'logFC'] < -1),]
up_phosp <- as.character(row.names(diff_up))
down_phosp <- as.character(row.names(diff_down))

pdf('/home/shengli/projects/AEG_proteomics/figures/Phosphoproteomics/diff_phosp_volcano_R1.pdf',width=6.87,height=5.92)
ggplot(diff,aes(x=logFC,y=lgfdr,label=phosp)) +
  geom_point(data=diff_none,aes(x=logFC,y=lgfdr),color='lightgray') +
  geom_point(data=diff_up,aes(x=logFC,y=lgfdr),color='red') +
  geom_point(data=diff_down,aes(x=logFC,y=lgfdr),color='blue') +
  geom_label_repel(data = subset(diff, lgfdr > 30 | abs(logFC) > 9),
                   size = 3,
                   box.padding = 0.5,
                   point.padding = 0.5,
                   segment.size = 0.2,
                   segment.color = "grey50") +
  theme(panel.background=element_rect(colour="black",fill="white"))
dev.off()

# barplot of functional enrichment analysis results, Fig. 7b
setwd('/home/shengli/projects/AEG_proteomics/results/phosphoproteome')
res_enrich <- read.table('GO_BP_Phosp_S3_diff.txt',header=T,sep='\t')
enrich_4bar <- data.frame(bp=res_enrich[1:8,'Description'],
                          fdr=res_enrich[1:8,'p.adjust'])
enrich_4bar$lgp <- -log10(enrich_4bar$fdr)
rownames(enrich_4bar) <- enrich_4bar[,'bp']
pdf('/home/shengli/projects/AEG_proteomics/figures/Phosphoproteomics/S3_diff_phosp_BP.pdf',width=7,height = 5)
barplot(enrich_4bar[,'lgp'],horiz=T,names.arg=rownames(enrich_4bar),las=2,xlim=c(0,16))
dev.off()

