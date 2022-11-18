### calculate difference of immune infiltration between paried tumor and normal samples

setwd('/home/shengli/projects/AEG_proteomics/results/RNAseq')
immune_infil <- read.table('Immune_CIBERSORT.txt',header=T,row.names=1,sep='\t')
immune_infil <- immune_infil[,1:22]

immune_samples <- rownames(immune_infil)
# S-I
samples_t1 <- c('C551594','C559459','C560016','C563560','C564447','C569091','C575352','C579894','C582702','C584162','C586163','C583851','C590054','C595498',
                'C596146','C1005614','C1006005','C1009828','C532443','C505919','C562401','C587488','C599521','C584104','C498355','C565023','C246567',
                'C576092','C488963','C533179','C504861','C568971','C528365','C533153','C1007477','C572176','C594667','C522528','C561723','C584243')
samples_t1_p <- samples_t1[which(samples_t1 %in% immune_samples)]
samples_n1 <- c('N551594','N559459','N560016','N563560','N564447','N569091','N575352','N579894','N582702','N584162','N586163','N583851','N590054','N595498',
                'N596146','N1005614','N1006005','N1009828','N532443','N505919','N562401','N587488','N599521','N584104','N498355','N565023','N246567',
                'N576092','N488963','N533179','N504861','N568971','N528365','N533153','N1007477','N572176','N594667','N522528','N561723','N584243')
samples_n1_p <- samples_n1[which(samples_n1 %in% immune_samples)]

Results <- matrix(NA,nrow=ncol(immune_infil),ncol=4,byrow=T)
rownames(Results) <- colnames(immune_infil)
colnames(Results) <- c('Tumor_infil','Normal_infil','infil_dif','Pvalue')

for (c in 1:ncol(immune_infil)) {
  infil_normal <- as.numeric(immune_infil[samples_n1_p,c])
  infil_tumor <- as.numeric(immune_infil[samples_t1_p,c])
  test_infil <- t.test(infil_normal,infil_tumor,paried=T)
  infil_ave_normal <- mean(infil_normal)
  infil_ave_tumor <- mean(infil_tumor)
  dif <- infil_ave_tumor - infil_ave_normal
  Results[c,'Tumor_infil'] <- infil_ave_tumor
  Results[c,'Normal_infil'] <- infil_ave_normal
  Results[c,'infil_dif'] <- dif
  Results[c,'Pvalue'] <- test_infil$p.value
}
tmp_pvalues <- as.numeric(Results[,'Pvalue'])
fdrs <- p.adjust(tmp_pvalues,method='fdr')
Results <- cbind(Immune_cell=colnames(immune_infil),Results,FDR=fdrs)
write.table(Results,file='Immune_infil_diff_g1.txt',row.names=F,quote=F,sep='\t')

# S-II
samples_t2 <- c('C217360','C266198','C325352','C227213','C502774','C284805','C538376','C335689','C242946','C242956','C505805','C194922','C197023','C195615',
                'C225468','C201779','C260052','C265603','C503345','C502439','C540265','C234760','C219685')
samples_t2_p <- samples_t2[which(samples_t2 %in% immune_samples)]
samples_n2 <- c('N217360','N266198','N325352','N227213','N502774','N284805','N538376','N335689','N242946','N242956','N505805','N194922','N197023','N195615',
                'N225468','N201779','N260052','N265603','N503345','N502439','N540265','N234760','N219685')
samples_n2_p <- samples_n2[which(samples_n2 %in% immune_samples)]
Results <- matrix(NA,nrow=ncol(immune_infil),ncol=4,byrow=T)
rownames(Results) <- colnames(immune_infil)
colnames(Results) <- c('Tumor_infil','Normal_infil','infil_dif','Pvalue')

for (c in 1:ncol(immune_infil)) {
  infil_normal <- as.numeric(immune_infil[samples_n2_p,c])
  infil_tumor <- as.numeric(immune_infil[samples_t2_p,c])
  test_infil <- t.test(infil_normal,infil_tumor,paried=T)
  infil_ave_normal <- mean(infil_normal)
  infil_ave_tumor <- mean(infil_tumor)
  dif <- infil_ave_tumor - infil_ave_normal
  Results[c,'Tumor_infil'] <- infil_ave_tumor
  Results[c,'Normal_infil'] <- infil_ave_normal
  Results[c,'infil_dif'] <- dif
  Results[c,'Pvalue'] <- test_infil$p.value
}
tmp_pvalues <- as.numeric(Results[,'Pvalue'])
fdrs <- p.adjust(tmp_pvalues,method='fdr')
Results <- cbind(Immune_cell=colnames(immune_infil),Results,FDR=fdrs)
write.table(Results,file='Immune_infil_diff_g2.txt',row.names=F,quote=F,sep='\t')

# S-III
samples_t3 <- c('C332654','C386485','C266473','C338321','C338421','C335207','C378018','C325830','C196713','C341159','C151058','C203659','C203933','C269944',
                'C210038','C260824','C218213','C391965','C363245','C196436','C321599','C241209','C215502','C248148','C322712','C248917','C0600290','C291416',
                'C1004644','C305692','C258083','C259302','C1003009','C288415','C314628','C260319','C199181','C256818','C260513','C260742')
samples_t3_p <- samples_t3[which(samples_t3 %in% immune_samples)]
samples_n3 <- c('N332654','N386485','N266473','N338321','N338421','N335207','N378018','N325830','N196713','N341159','N151058','N203659','N203933','N269944',
                'N210038','N260824','N218213','N391965','N363245','N196436','N321599','N241209','N215502','N248148','N322712','N248917','N0600290','N291416',
                'N1004644','N305692','N258083','N259302','N1003009','N288415','N314628','N260319','N199181','N256818','N260513','N260742')
samples_n3_p <- samples_n3[which(samples_n3 %in% immune_samples)]
Results <- matrix(NA,nrow=ncol(immune_infil),ncol=4,byrow=T)
rownames(Results) <- colnames(immune_infil)
colnames(Results) <- c('Tumor_infil','Normal_infil','infil_dif','Pvalue')

for (c in 1:ncol(immune_infil)) {
  infil_normal <- as.numeric(immune_infil[samples_n3_p,c])
  infil_tumor <- as.numeric(immune_infil[samples_t3_p,c])
  test_infil <- t.test(infil_normal,infil_tumor,paried=T)
  infil_ave_normal <- mean(infil_normal)
  infil_ave_tumor <- mean(infil_tumor)
  dif <- infil_ave_tumor - infil_ave_normal
  Results[c,'Tumor_infil'] <- infil_ave_tumor
  Results[c,'Normal_infil'] <- infil_ave_normal
  Results[c,'infil_dif'] <- dif
  Results[c,'Pvalue'] <- test_infil$p.value
}
tmp_pvalues <- as.numeric(Results[,'Pvalue'])
fdrs <- p.adjust(tmp_pvalues,method='fdr')
Results <- cbind(Immune_cell=colnames(immune_infil),Results,FDR=fdrs)
write.table(Results,file='Immune_infil_diff_g3.txt',row.names=F,quote=F,sep='\t')

