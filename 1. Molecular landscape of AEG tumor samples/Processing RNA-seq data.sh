#!/bin/bash
#hisat2 version 2.2.1
#reference human genome GENCODE v35
#sambamba version 0.8.0
#samtools version 1.10

sample=$1
dir_wk='/extraspace/projects/AEG_proteomics'
dir_fq=${dir_wk}'/fastq'
fq1=${dir_fq}/${sample}/${sample}'_R1.fastq.gz'
fq2=${dir_fq}/${sample}/${sample}'_R2.fastq.gz'
ht2='/home/Softwares/hisat2-2.2.1/hisat2'
ht2_indx='/home/Reference/fa/human/hisat2_gencodev35_grch38/hisat2_gencodev35_grch38'
outdir=${dir_wk}'/alignmet_ht2/'${sample}
mkdir -p ${outdir}
out_sam=${outdir}/${sample}'_alignment.sam'
${ht2} -x ${ht2_indx} -1 ${fq1} -2 ${fq2} -S ${out_sam}
sambamba='/home/Softwares/sambamba-0.8.0-linux-amd64-static'
bamfile=${outdir}/${sample}'_alignment.bam'
${sambamba} view ${out_sam} -S -f bam -o ${bamfile}
samtools='/usr/bin/samtools'
file_stat=${dir_wk}'/alignment_ht2/'${sample}'/'${sample}'_alignment.stat'
${samtools} flagstat ${bamfile} > ${file_stat}
sortbam=${dir_wk}'/alignment_ht2/'${sample}'/'${sample}'_alignment_sorted.bam'
tmpdir=${dir_wk}'/alignment_ht2/'${sample}
${sambamba} sort -o ${sortbam} --tmpdir ${tmpdir} -t 8 ${bamfile}
rm -f ${out_sam}
rm -f ${bamfile}

