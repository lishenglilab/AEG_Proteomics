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

## somatic interactions analysis
rm(list=ls())
library(maftools)
setwd('/home/shengli/projects/AEG_proteomics/results/WES')

# S-I subtype
aeg1_mut <- read.maf('AEG_mutation_neat_filter_non151058_group1.maf')
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/SomaticInteractions_g1.pdf',width=5, height=5)
somaticInteractions(maf=aeg1_mut,top=25,pvalue=c(0.05,0.1))
dev.off()

# S-II subtype
aeg2_mut <- read.maf('AEG_mutation_neat_filter_non151058_group2.maf')
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/SomaticInteractions_g2.pdf',width=5, height=5)
somaticInteractions(maf=aeg2_mut,top=25,pvalue=c(0.05,0.1))
dev.off()

# S-III subtype
aeg3_mut <- read.maf('AEG_mutation_neat_filter_non151058_group3.maf')
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/SomaticInteractions_g3.pdf',width=5, height=5)
somaticInteractions(maf=aeg3_mut,top=25,pvalue=c(0.05,0.1))
dev.off()

## Oncogenic pathway analysis
rm(list=ls())
library(maftools)
library(R.utils)
setwd('/home/shengli/projects/AEG_proteomics/results/WES')

#group1
aeg1_mut <- read.maf('AEG_mutation_neat_filter_non151058_group1.maf')
OncogenicPathways(maf = aeg1_mut)
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/OncogenicPathway/NOTCH_pathway_g1.pdf',width=6,height = 5)
PlotOncogenicPathways(maf = aeg1_mut, pathways = "NOTCH")
dev.off()
#group2
aeg2_mut <- read.maf('AEG_mutation_neat_filter_non151058_group2.maf')
OncogenicPathways(maf = aeg2_mut)
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/OncogenicPathway/RTK-RAS_pathway_g2.pdf',width=6,height = 5)
PlotOncogenicPathways(maf = aeg2_mut, pathways = "RTK-RAS")
dev.off()
#group3
aeg3_mut <- read.maf('AEG_mutation_neat_filter_non151058_group3.maf')
OncogenicPathways(maf = aeg3_mut)
pdf('/home/shengli/projects/AEG_proteomics/figures/WES/OncogenicPathway/RTK-RAS_pathway_g3.pdf',width=6,height = 5)
PlotOncogenicPathways(maf = aeg3_mut, pathways = "RTK-RAS")
dev.off()

### calculate the effect of mutations on proteome and phosphoproteome

setwd('/home/shengli/projects/AEG_proteomics/results/proteome')
prot_fc <- read.table('Proteomics_iBAQ103_log2quantile_normlization_impute_perc25_TN.txt',header=T,as.is=1,row.names=1)

setwd('/home/shengli/projects/AEG_proteomics/results/WES/2022.08')

mut_samples <- read.table('AEG_mut_samples.txt',header=T,sep='\t')
genes <- c('LEPR','SORCS3','POTEE','COL2A1','DLG2','PHIP','ANKRD36C','MAG','TBX22','RPS6KA6','RP11-267C16.1','PAX8','DIDO1','TRPC6','KCNH2','PDGFC','STAB2',
           'ZNF831','NOTCH3','ADAM22','PTPRQ','NRK','TP53','DNM1P47','FAT3','REG3A','MAD1L1','SETD1B','SNURF','ACSM2A','KAT2A','LILRB3','XRN2','USP4','NEDD9',
           'KLF12','TRIP11','ATXN2L','ASMTL','TNRC18P2','CCNA1','RLTPR','RP11-726G1.1','DHX35','DSC1','PTPN23','TMEM225','A2ML1','ANTXR1','ZNF385B','ITGA2B',
           'SGCD','SEPT12','MROH2B','PRKG1','SAMD3','C18orf64','BPI','GUCY1A2','HERC2','PREX2','AC024937.6','COL24A1','SLIT2','SPTA1','MUC19','MRS2','SDK1',
           'COL14A1','MYOF','RRP12','HIP1R','HELZ2','MTHFD1L','INSC','USP7','NCOR1','AXDND1','ADCY5','ENPEP','ELMO1','RGS6','CYFIP1','NOS3','CACNA1F','COL4A4',
           'COL19A1','ALPI','COL13A1','GFRA1','ESRRG','VPS13C','ZFHX4','ABI3BP','DPP10','GRID2','AMPH','NCKAP1','DOPEY1','DENND1B','TRPC4','MAP7','CPED1','HNRNPR',
           'NR1H3','GREB1','PALD1','LHX9','PPP5C','TRIM6','SLITRK1','ARID2','SMAD4','AC024937.4','POSTN','RPL23AP7','GALNS','CDHR2','COLGALT1','FAM83H',
           'CACNA1S','SLC12A6','SON','DTNA','NUP98','OVGP1','ATP5L','SNAI2','SUPV3L1','COL5A1','ZNF99','WDR26','DAAM2','MACF1','NAV2','RP11-166O4.1','MST1L',
           'RP11-782C8.2','C9orf172','MYBPC3','RAD9A','ZFC3H1','CIITA','DOK5','NCAPH2','ALMS1','SLC9B1','SUN1','NLRP11','NGFRAP1','GORASP1','ZNF595','LRWD1',
           'MTRF1','UGGT2','COG7','MFI2','PRSS12','NAA16','FCGR3A','VPS54','ANKRD42','KREMEN2','MIS18A','MAPK14','KEL','DAAM1','THUMPD1','MARK4','TXLNG','CAMSAP3',
           'DACT3','GART','ORC2','STK38L','C1orf109','KIAA0754','BCO2','TUBB3','RIPK1','EPM2A','DPP7','OR4C5','RIN1','RDX','KLRB1','SDSL','MARK1','OSMR','ARHGAP8',
           'LACTBL1','KCNQ5','CLN6','B4GALNT1','EPS8L3','AIDA','ZNF721','NANOGP1','TUBG2','CCDC8','DSCC1','SAT2','FBXO31','SLC35E3','CCHCR1','MUC15','RB1','SLC9A6',
           'SUPT3H','NAT14','RPLP1','RBM41','ZBTB5','RASSF8','RAB11FIP3','PTX4','SPN','PTP4A3','OR1D4','ZNF691','EXOSC4','ZDHHC21','DNAJC5G','U2SURP','DTX3','RASD2',
           'TTTY14','VARS2','CHMP2B','RNF112','ZDHHC12','TRIM51HP','CPO','CPB1','OR4X1','MAGED2','TRBV2','TAS2R13','IGHV3-74','PPME1','WIZ','LGALS9','CALCRL','FGF18',
           'REEP6','WDR46','HMGXB3','ZZEF1','ST3GAL3','USP9X','PTPRB','TMTC4','ADRA1D','LAS1L','ADAMTS3','CSMD2','ANKRD20A9P','BCORL1','THSD7A','TAF1','SLC51A','SPTBN2',
           'TLR9','CLCNKA','USH2A','CTD-2089O24.1','PXDN','DST','SVEP1'
           )

samples_prot <- colnames(prot_fc)

gp_pairs <- c()
pvalues <- c()
changes <- c()
for (g in genes) {
  samples_g <- as.character(mut_samples[which(mut_samples[,'Gene'] == g),'Sample'])
  samples_g <- unique(samples_g)
  samples_ng <- samples_prot[!(samples_prot %in% samples_g)]
  for (prt in rownames(prot_fc)) {
    expr_mut <- as.numeric(prot_fc[prt,samples_g])
    expr_wt <- as.numeric(prot_fc[prt,samples_ng])
    change_mut_wt <- mean(expr_mut) - mean(expr_wt)
    comp <- t.test(log2(expr_mut),log2(expr_wt))
    if (comp$p.value == 'NaN') {
      next
    }
    pair <- paste(g,prt,sep='\t')
    gp_pairs <- c(gp_pairs,pair)
    pvalues <- c(pvalues,comp$p.value)
    changes <- c(changes,change_mut_wt)
  }
}
fdrs <- p.adjust(pvalues,method='fdr')
res <- cbind(gp_pairs,changes,pvalues,fdrs)
colnames(res) <- c('Gene  Protein','FoldChange','Pvalue','FDR')
setwd('/home/shengli/projects/AEG_proteomics/results/WES/2022.08')
write.table(res,file='AEG_specific_mut_affect_proteins.txt',sep='\t',quote=F,row.names=F)

