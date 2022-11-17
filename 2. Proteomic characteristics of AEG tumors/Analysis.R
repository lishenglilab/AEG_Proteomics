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

###

