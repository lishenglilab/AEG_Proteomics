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

