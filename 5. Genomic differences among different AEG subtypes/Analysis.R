### mutation signature analysis
rm(list=ls())
library(maftools)
library(BSgenome)
library('NMF')
library(BSgenome.Hsapiens.UCSC.hg19)
setwd('/home/shengli/projects/AEG_proteomics/results/WES')
#S-I subtype
aeg1_mut <- read.maf('AEG_mutation_neat_filter_non151058_group1.maf')
aeg1_tnm <- trinucleotideMatrix(aeg1_mut,ref_genome='BSgenome.Hsapiens.UCSC.hg19')
aeg1_sign <- estimateSignatures(mat=aeg1_tnm,nTry=6)
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/estimate_mutSig_cophenetic_g1.pdf',width=5,height=5)
plotCophenetic(res=aeg1_sign)
dev.off()

aeg1_sig <- extractSignatures(mat = aeg1_tnm, n = 3)
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/mutSig_SBS_g1.pdf',width=5,height = 6)
plotSignatures(nmfRes = aeg1_sig, sig_db = 'SBS')
dev.off()

# S-II subtype
aeg2_mut <- read.maf('AEG_mutation_neat_filter_non151058_group2.maf')
aeg2_tnm <- trinucleotideMatrix(aeg2_mut,ref_genome='BSgenome.Hsapiens.UCSC.hg19')
aeg2_sign <- estimateSignatures(mat=aeg2_tnm,nTry=6)
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/estimate_mutSig_cophenetic_g2.pdf',width=5,height=5)
plotCophenetic(res=aeg2_sign)
dev.off()

aeg2_sig <- extractSignatures(mat = aeg2_tnm, n = 3)
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/mutSig_SBS_g2.pdf',width=5,height = 6)
plotSignatures(nmfRes = aeg2_sig, sig_db = 'SBS')
dev.off()
# S-III subtype
aeg3_mut <- read.maf('AEG_mutation_neat_filter_non151058_group3.maf')
aeg3_tnm <- trinucleotideMatrix(aeg3_mut,ref_genome='BSgenome.Hsapiens.UCSC.hg19')
aeg3_sign <- estimateSignatures(mat=aeg3_tnm,nTry=6)
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/estimate_mutSig_cophenetic_g3.pdf',width=5,height=5)
plotCophenetic(res=aeg3_sign)
dev.off()

aeg3_sig <- extractSignatures(mat = aeg3_tnm, n = 3)
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/mutSig_SBS_g3.pdf',width=5,height = 6)
plotSignatures(nmfRes = aeg3_sig,sig_db = 'SBS')
dev.off()
