### Volcano plot showing the difference in proteins between AEG tumor and paired NAT samples. Fig. 2a
## calculate the difference
library(limma) # version 3.46.0
setwd('/home/shengli/projects/0.Collaborations/4.Qin.Jiangjiang/data/proteome')
prot_mx <- read.table('Proteomics_iBAQ103_log2quantile_normlization_impute_perc25.txt',header=T,row.names=1,sep='\t')
design_nt <- cbind(Intercept=1,
                   Group=c(rep(1,103),rep(0,103)))
fit <- lmFit(prot_mx,design_nt)
fit <- eBayes(fit,trend=TRUE)
dep_nt <- topTable(fit,coef=2,number=Inf)
setwd('/home/shengli/projects/0.Collaborations/4.Qin.Jiangjiang/results/proteome')
write.table(dep_nt,file="Normal_Tumor_DEP_all.txt",quote=F,sep='\t',append=F)

### Functional enrichment results of up-regulated and down-regulated proteins, respectively. Fig. 2b

### Fig. 2c
