#!/bin/bash

### estimate the optimal rank

r_estimate='/home/Softwares/RBPgroup/bin/NMF.estimate.R'
file_mx='/home/shengli/projects/AEG_proteomics/results/Proteomics_perc25_impute_cvs25.txt'
dir_est='/home/shengli/projects/AEG_proteomics/results/NMF/output/estiRank'
${r_estimate} -i ${file_mx} -s 2 -e 7 -n 200 -o ${dir_est}/Proteomics.estiRank

### Run NMF with the selected rank
r_main='/home/Softwares/RBPgroup/bin/NMF.main.R'
file_mx='/home/shengli/projects/AEG_proteomics/results/Proteomics_perc25_impute_cvs25.txt'
dir_main='/home/shengli/projects/AEG_proteomics/results/NMF/output/main'
${r_main} -i ${file_mx} -r 3 -n 200 -o ${dir_main}/Proteomics.main3

### Extract the cluster components

r_assign='/home/Softwares/RBPgroup/bin/NMF.assign_clusters.R'
file_coef='/home/shengli/projects/AEG_proteomics/results/NMF/output/main/Proteomics.main2.coef'
dir_assign='/home/shengli/projects/AEG_proteomics/results/NMF/output/assign_clusters'
${r_assign} -i ${file_coef} -o ${dir_assign}/Proteomics.2.assign_cluster.txt
