### waterfall plot of top mutation cancer genes, Fig. 1b
library(maftools)
setwd('/home/shengli/projects/AEG_proteomics/data/WES')
aeg <- read.maf(maf='AEG_mutation_cancer_gene.maf')
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/mutation_waterfall_cancer_genes_top30.pdf',height=6,width=10)
oncoplot(maf=aeg,draw_titv=T,top=30,showTumorSampleBarcodes=F)
dev.off()

###
