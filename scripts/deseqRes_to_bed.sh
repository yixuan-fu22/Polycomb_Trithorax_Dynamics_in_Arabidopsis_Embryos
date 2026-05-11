#!/bin/bash

awk 'BEGIN{OFS="\t"} NR>1 {print $2, $3, $4, $12}' DESeq2Res_withAnnotations.tsv > DEseqRes.bed

# $12 = FDR

# batch, resutls folder

for i in $(ls -d diffBind_ChIPseeker*); do 
awk 'BEGIN{OFS="\t"} NR>1 {print $2, $3, $4, $12}' ${i}/DESeq2Res_withAnnotations.tsv > ${i}/${i/diffBind_ChIPseeker_/}_DEseq2Res.bed; done

#

awk 'BEGIN{OFS="\t"} $2 != "ChrC" && $2 != "ChrM" && NR>1 && $12 < 0.05 && $8 > $9 {print $2,$3,$4,$12}' DESeq2Res_withAnnotations.tsv > diffPeaks_bed_upInStage1.bed

awk 'BEGIN{OFS="\t"} $2 != "ChrC" && $2 != "ChrM" && NR>1 && $12 < 0.05 && $8 < $9 {print $2,$3,$4,$12}' DESeq2Res_withAnnotations.tsv > diffPeaks_bed_upInStage2.

#

for i in $(ls -d diffBind_ChIPseeker*); do 
awk 'BEGIN{OFS="\t"} $2 != "ChrC" && $2 != "ChrM" && NR>1 && $12 < 0.05 && $8 > $9 {print $2,$3,$4,$12}' ${i}/DESeq2Res_withAnnotations.tsv > ${i}/${i/"diffBind_ChIPseeker_"/}_diffPeaks_upInStage1.bed; 
awk 'BEGIN{OFS="\t"} $2 != "ChrC" && $2 != "ChrM" && NR>1 && $12 < 0.05 && $8 < $9 {print $2,$3,$4,$12}' ${i}/DESeq2Res_withAnnotations.tsv > ${i}/${i/"diffBind_ChIPseeker_"/}_diffPeaks_upInStage2.bed
done