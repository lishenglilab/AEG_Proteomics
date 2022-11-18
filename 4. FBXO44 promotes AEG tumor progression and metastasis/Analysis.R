### calculate the associations between FBXO44 and immune checkpoint genes
setwd('/home/shengli/projects/AEG_proteomics/results/RNAseq/Immune')
gene_expr <- read.table('AEG_immck_gene_profile.txt',header=T,row.names=1,sep='\t')
gene_expr <- gene_expr[,104:206]
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
### calculate the associations between FBXO44 and immune cell abundance
setwd('/home/shengli/projects/AEG_proteomics/results/RNAseq/Immune')
gene_expr <- read.table('AEG_immck_gene_profile.txt',header=T,row.names=1,sep='\t')
gene_expr <- gene_expr[,104:206]
fbxo44_expr <- as.numeric(gene_expr['FBXO44',])
names(fbxo44_expr) <- colnames(gene_expr)
cell_abd <- read.table('xCell_AEG_gene_tmp_mx_ordered_xCell_2052121421.txt',header=T,row.names=1,sep='\t')
cells <- c('Plasma cells','CD4+ Tem','Th1 cells','Mast cells','CD4+ memory T-cells','Tregs','Macrophages M2','naive B-cells','cDC','Macrophages',
           'Basophils','Macrophages M1','Tgd cells','iDC','DC','aDC','Monocytes','CD4+ Tcm','Th2 cells','pDC')
cell_abd <- cell_abd[cells,104:206]

fbxo44_expr_ord <- fbxo44_expr[colnames(cell_abd)]
res_pvalues <- c()
res_cors <- c()
for (n in 1:nrow(cell_abd)) {
  abd <- as.numeric(cell_abd[n,])
  cor_test <- cor.test(fbxo44_expr_ord,abd)
  res_pvalues <- c(res_pvalues,cor_test$p.value)
  res_cors <- c(res_cors,cor_test$estimate)
}
res_fdrs <- p.adjust(res_pvalues,method='fdr')
res <- cbind(rownames(cell_abd),res_cors,res_pvalues,res_fdrs)
colnames(res) <- c('Cell','Correlation','Pvalue','FDR')
write.table(res, file='FBXO44_immunecell_cor.txt',quote=F,sep='\t',row.names=F)


