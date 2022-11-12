### waterfall plot of top mutation cancer genes, Fig. 1b
library(maftools)
setwd('/home/shengli/projects/AEG_proteomics/data/WES')
aeg <- read.maf(maf='AEG_mutation_cancer_gene.maf')
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/mutation_waterfall_cancer_genes_top30.pdf',height=6,width=10)
oncoplot(maf=aeg,draw_titv=T,top=30,showTumorSampleBarcodes=F)
dev.off()

### The number distribution of identified proteins in each sample, Fig. 1c
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
prt_num <- read.table('Proteomics_sample_numbers.txt',header=T,sep='\t')
num_tumor <- prt_num[which(prt_num[,'Type']=='tumor'),]
rownames(num_tumor) <- num_tumor[,'Sample']
num_tumor_ordered <- num_tumor[order(num_tumor[,'Protein_number'],decreasing = F),]
pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics/Protein_numbers_tumor.pdf',height=5,width=5)
plot(num_tumor_ordered[,'Protein_number'],ylim=c(8000,9500),las=1)
dev.off()
num_normal <- prt_num[which(prt_num[,'Type']=='normal'),]
rownames(num_normal) <- num_normal[,'Sample']
num_normal_ordered <- num_normal[rownames(num_tumor_ordered),]
wilcox.test(as.numeric(num_normal_ordered[,'Protein_number']),as.numeric(num_tumor_ordered[,'Protein_number']))

pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics/Protein_numbers_normal.pdf',height=5,width=5)
plot(num_normal_ordered[,'Protein_number'],ylim=c(8000,9500),las=1)
dev.off()
library(ggplot2)
nums <- c(as.numeric(num_tumor_ordered[,'Protein_number']),as.numeric(num_normal_ordered[,'Protein_number']))
nums_df <- data.frame(Numbers=nums,
                      sample_group=c(rep('tumor',103),rep('normal',103)))
pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics/protein_numbers_comp.pdf',height=6,width=5)
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

### The number distribution of phosphorylation sites in each sample, Fig. 1d
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

### The number distribution of expressed genes in each sample, Fig. 1e
setwd('/home/shengli/projects/AEG_proteomics/results/RNAseq')
gene_num <- read.table('gene_sample_numbers.txt',header=T,sep='\t')

num_normal <- gene_num[which(gene_num[,'Type']=='normal'),]
rownames(num_normal) <- num_normal[,'Sample']
num_normal_ordered <- num_normal[order(num_normal[,'Gene_number'],decreasing = F),]
pdf('/home/shengli/projects/AEG_proteomics/figures/RNAseq/Gene_numbers_normal.pdf',height=5,width=5)
barplot(num_normal_ordered[,'Gene_number'],ylim=c(0,30000),las=1)
dev.off()

num_tumor <- gene_num[which(gene_num[,'Type']=='tumor'),]
rownames(num_tumor) <- num_tumor[,'Sample']
num_tumor_ordered <- num_tumor[rownames(num_normal_ordered),]
pdf('/home/shengli/projects/AEG_proteomics/figures/RNAseq/Gene_numbers_tumor.pdf',height=5,width=5)
barplot(num_tumor_ordered[,'Gene_number'],ylim=c(0,30000),las=1)
dev.off()


### Boxplot showing the comparison of the numbers of expressed genes between parired AEG tumor and NAT samples, Supplementary Fig. 4c
library(ggplot2)
nums <- c(as.numeric(num_tumor_ordered[,'Gene_number']),as.numeric(num_normal_ordered[,'Gene_number']))
nums_df <- data.frame(Numbers=nums,
                      sample_group=c(rep('tumor',85),rep('normal',85)))
pdf('/home/shengli/projects/AEG_proteomics/figures/RNAseq/gene_numbers_comp.pdf',height=6,width=5)
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
wilcox.test(as.numeric(nums_df[which(nums_df[,'sample_group']=='normal'),'Numbers']),as.numeric(nums_df[which(nums_df[,'sample_group']=='tumor'),'Numbers']))


