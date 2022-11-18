### calculate the associations between FBXO44 and immune checkpoint genes
setwd('/home/shengli/projects/AEG_proteomics/results/RNAseq/Immune')
gene_expr <- read.table('AEG_immck_gene_profile.txt',header=T,row.names=1,sep='\t')
gene_expr <- gene_expr[,104:206]
# immune checkpoints
fbxo44_expr <- as.numeric(gene_expr['FBXO44',])
res_pvalues <- c()
res_cors <- c()
for (n in 1:nrow(gene_expr)) {
  expr <- as.numeric(gene_expr[n,])
  cor_test <- cor.test(expr,fbxo44_expr)
  res_pvalues <- c(res_pvalues,cor_test$p.value)
  res_cors <- c(res_cors,cor_test$estimate)
}
res_fdrs <- p.adjust(res_pvalues,method='fdr')
res <- cbind(rownames(gene_expr),res_cors,res_pvalues,res_fdrs)
colnames(res) <- c('Gene','Correlation','Pvalue','FDR')
write.table(res,file='FBXO44_immunecheck_cor.txt',quote=F,sep='\t',row.names=F)


