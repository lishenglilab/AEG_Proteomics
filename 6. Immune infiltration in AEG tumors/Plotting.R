### heatmap of immune cells in three AEG subtypes, Fig. 6a
setwd('/home/shengli/projects/AEG_proteomics/results/RNAseq/Immune')
cell_abd <- read.table('AEG_xCell_subtypes_zscore.txt',header=T,sep='\t')

cells_ord <- c('Adipocytes','Astrocytes','Fibroblasts','Preadipocytes','Pericytes','ly Endothelial cells','mv Endothelial cells',
               'Smooth muscle','Endothelial cells','Myocytes','Basophils','DC','Macrophages','Macrophages M1',
               'Macrophages M2','Mast cells','aDC','cDC','iDC','pDC','Monocytes','Platelets','CLP','CMP','GMP','HSC','MEP','Megakaryocytes',
               'Epithelial cells','Mesangial cells','Neurons','Sebocytes',
               'Tregs','Tgd cells','CD4+ Tcm','CD4+ Tem','CD4+ memory T-cells','Th1 cells','Th2 cells','naive B-cells','Plasma cells') # lymphoid cells(9), myeloid cells(11), stromal cells(10), stem cells(7), other cells(4), 
samples_ord <- c(samples_t1,samples_t2,samples_t3)
pdf('/home/shengli/projects/AEG_proteomics/figures/RNAseq/xCell_subtype_heatmap.pdf',height=8,width=8)
ggplot(cell_abd,aes(x=Sample,y=Cell)) +
  geom_tile(aes(fill=Abundance),col='lightgray') +
  scale_fill_gradient2(low='blue',high='red')
  scale_x_discrete(limit=samples_ord) +
  scale_y_discrete(limit=cells_ord) +
  theme(panel.background = element_rect(colour=NA,fill='white',size=2),
        panel.grid = element_blank(),
        panel.grid.major=element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_text(size=10,colour='black',angle=90,hjust=1,vjust=0.5),
        axis.text.y=element_text(size=10,colour='black'),
        axis.ticks=element_blank())
dev.off()

# difference between subtypes
sub_pvalues <- c()
for (n in 1:nrow(cell_abd)) {
  abd_s1 <- as.numeric(cell_abd[n,samples_t1])
  abd_s2 <- as.numeric(cell_abd[n,samples_t2])
  abd_s3 <- as.numeric(cell_abd[n,samples_t3])
  abd_sub_df <- data.frame(Abundance = c(abd_s1, abd_s2, abd_s3),
                           Group = c(rep(1,length(abd_s1)),rep(2,length(abd_s2)),rep(3,length(abd_s3))))
  comp <- kruskal.test(Abundance~Group,data=abd_sub_df)
  sub_pvalues <- c(sub_pvalues,comp$p.value)
}
res_sub <- cbind(rownames(cell_abd),sub_pvalues)

# plot individual cells, Supplementary Fig. 14
cell_abd_tumor <- cell_abd[,samples_c1c2c3]

cel_abd <- as.numeric(cell_abd_tumor['mv Endothelial cells',c(samples_c1c2c3)])

cel_abd_df <- data.frame(Abundance = cel_abd,
                         Group = c(rep('S1',length(samples_t1)),rep('S2',length(samples_t2)),rep('S3',length(samples_t3))))

pdf('/home/shengli/projects/AEG_proteomics/figures/RNAseq/xCell_sub_box_mv Endothelial cells.pdf',height = 6,width = 3)
ggplot(cel_abd_df,aes(x=Group,y=Abundance)) +
  geom_jitter(width=0.15,size=0.7) +
  geom_boxplot(color="black",fill=NA,outlier.shape=NA,width=0.3) +
  scale_x_discrete(limit=c('S1','S2','S3'),expand=c(0.1,0)) +
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

### plot differences of immune cells, Fig. 6b
library(ggplot2)
setwd('/home/shengli/projects/AEG_proteomics/results/RNAseq')
immune_dif <- read.table('Immune_infil_xCell_difference_df.txt',header=T,sep='\t')
immune_dif$lgp <- -log10(immune_dif[,'Pvalue'])
immune_cells <- immune_dif[,'Immune_cell']
pdf('/home/shengli/projects/AEG_proteomics/figures/RNAseq/Immune_infil_xCell_diff.pdf',width=9.5,height=3.5)
ggplot(immune_dif,aes(x=Immune_cell,y=Group))+
  geom_point(aes(size=lgp,col=Difference))+
  scale_color_gradient2(low="blue",high="red",na.value="white",name="Difference")+
  scale_size_continuous(limit=c(0,15),range=c(0.5,4),breaks=c(-log10(0.05),5,10),labels=c("0.05","1e-5","1e-10"))+
  scale_y_discrete(limit=c('Group1','Group2','Group3'),expand=c(0.05,0.05))+
  scale_x_discrete(limit=immune_cells,expand=c(0.05,0.05)) +
  theme(panel.background=element_rect(colour="black",fill="white"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=11,colour="black"),
        axis.text.x=element_text(size=10,colour="black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.key=element_rect(fill="white",colour="black"))
dev.off()

# plot individual immune cells, Fig. 6c-d
library(ggplot2)
setwd('/home/shengli/projects/AEG_proteomics/results/RNAseq')
immune_prof <- read.table('xCell_AEG_gene_tmp_mx_ordered_xCell_2052121421.txt',header=T,row.names=1,sep='\t')
immune_samples <- colnames(immune_prof)
samples_t1 <- c('C551594','C559459','C560016','C563560','C564447','C569091','C575352','C579894','C582702','C584162','C586163','C583851','C590054','C595498',
                'C596146','C1005614','C1006005','C1009828','C532443','C505919','C562401','C587488','C599521','C584104','C498355','C565023','C246567',
                'C576092','C488963','C533179','C504861','C568971','C528365','C533153','C1007477','C572176','C594667','C522528','C561723','C584243')
samples_t1_p <- samples_t1[which(samples_t1 %in% immune_samples)]
samples_n1 <- c('N551594','N559459','N560016','N563560','N564447','N569091','N575352','N579894','N582702','N584162','N586163','N583851','N590054','N595498',
                'N596146','N1005614','N1006005','N1009828','N532443','N505919','N562401','N587488','N599521','N584104','N498355','N565023','N246567',
                'N576092','N488963','N533179','N504861','N568971','N528365','N533153','N1007477','N572176','N594667','N522528','N561723','N584243')
samples_n1_p <- samples_n1[which(samples_n1 %in% immune_samples)]
samples_t2 <- c('C217360','C266198','C325352','C227213','C502774','C284805','C538376','C335689','C242946','C242956','C505805','C194922','C197023','C195615',
                'C225468','C201779','C260052','C265603','C503345','C502439','C540265','C234760','C219685')
samples_t2_p <- samples_t2[which(samples_t2 %in% immune_samples)]
samples_n2 <- c('N217360','N266198','N325352','N227213','N502774','N284805','N538376','N335689','N242946','N242956','N505805','N194922','N197023','N195615',
                'N225468','N201779','N260052','N265603','N503345','N502439','N540265','N234760','N219685')
samples_n2_p <- samples_n2[which(samples_n2 %in% immune_samples)]
samples_t3 <- c('C332654','C386485','C266473','C338321','C338421','C335207','C378018','C325830','C196713','C341159','C151058','C203659','C203933','C269944',
                'C210038','C260824','C218213','C391965','C363245','C196436','C321599','C241209','C215502','C248148','C322712','C248917','C0600290','C291416',
                'C1004644','C305692','C258083','C259302','C1003009','C288415','C314628','C260319','C199181','C256818','C260513','C260742')
samples_t3_p <- samples_t3[which(samples_t3 %in% immune_samples)]
samples_n3 <- c('N332654','N386485','N266473','N338321','N338421','N335207','N378018','N325830','N196713','N341159','N151058','N203659','N203933','N269944',
                'N210038','N260824','N218213','N391965','N363245','N196436','N321599','N241209','N215502','N248148','N322712','N248917','N0600290','N291416',
                'N1004644','N305692','N258083','N259302','N1003009','N288415','N314628','N260319','N199181','N256818','N260513','N260742')
samples_n3_p <- samples_n3[which(samples_n3 %in% immune_samples)]
cell <- c('aDC')
s1_normal = as.numeric(immune_prof[cell,samples_n1_p])
s1_tumor = as.numeric(immune_prof[cell,samples_t1_p])
s2_normal = as.numeric(immune_prof[cell,samples_n2_p])
s2_tumor = as.numeric(immune_prof[cell,samples_t2_p])
s3_normal = as.numeric(immune_prof[cell,samples_n3_p])
s3_tumor = as.numeric(immune_prof[cell,samples_t3_p])
cell_abundance <- c(s1_normal,s1_tumor,s2_normal,s2_tumor,s3_normal,s3_tumor)
cell_abund <- data.frame(Abundance = cell_abundance,
                         Group = c(rep('S1_normal',length(s1_normal)),rep('S1_tumor',length(s1_tumor)),rep('S2_normal',length(s2_normal)),
                                   rep('S2_tumor',length(s2_tumor)),rep('S3_normal',length(s3_normal)),rep('S3_tumor',length(s3_tumor))))
pdf('/home/shengli/projects/AEG_proteomics/figures/RNAseq/xCell_cell_abundance_aDC.pdf',height=5,width=6)
ggplot(cell_abund,aes(x=Group,y=Abundance)) +
  geom_jitter(width=0.15,size=0.7) +
  geom_boxplot(color="black",fill=NA,outlier.shape=NA,width=0.3) +
  scale_x_discrete(limit=c('S1_normal','S1_tumor','S2_normal','S2_tumor','S3_normal','S3_tumor'),expand=c(0.1,0)) +
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
# plot chekpoint difference, Fig. 6f
library(ggplot2)
setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
ckp_diff <- read.table('G1G2G3_NAT_DEP_checkpoint_df.txt',header=T,sep='\t')
ckp_diff$lgp <- -log10(ckp_diff[,'FDR'])
ckp_diff[which(ckp_diff[,'lgp'] > 15),'lgp'] <- 15
ckp_diff[which(ckp_diff[,'logFC'] > 5),'logFC'] <- 5

pdf('/home/shengli/projects/AEG_proteomics/figures/Proteomics/Checkpoint_DEP_diff.pdf',width=9.5,height=3.5)
ggplot(ckp_diff,aes(x=Gene,y=Group))+
  geom_point(aes(size=lgp,col=logFC))+
  scale_color_gradient2(low="blue",high="red",na.value="white",name="Difference")+
  scale_size_continuous(limit=c(0,15),range=c(0.5,4),breaks=c(-log10(0.05),5,10),labels=c("0.05","1e-5","1e-10"))+
  scale_y_discrete(limit=c('G1','G2','G3'),expand=c(0.05,0.05))+
  theme(panel.background=element_rect(colour="black",fill="white"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=11,colour="black"),
        axis.text.x=element_text(size=10,colour="black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.key=element_rect(fill="white",colour="black"))
dev.off()


