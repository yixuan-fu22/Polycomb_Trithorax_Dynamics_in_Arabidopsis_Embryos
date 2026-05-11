#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate bedtools

peakdir="../data_fromOtherProjects/bedtools_consensusPeak"
outdir="../results/substractedPeaks"
mkdir -p ${outdir}

# H2Aub
bedtools subtract -A -a ${peakdir}/WS_H2AUB.consensusPeaks.bed \
-b ${peakdir}/MG_H2AUB.consensusPeaks.bed > ${outdir}/H2AUB_WS_substract_MG.bed

bedtools subtract -A -a ${peakdir}/MG_H2AUB.consensusPeaks.bed \
-b ${peakdir}/WS_H2AUB.consensusPeaks.bed > ${outdir}/H2AUB_MG_substract_WS.bed

cat ${peakdir}/WS_H2AUB.consensusPeaks.bed \
 ${peakdir}/MG_H2AUB.consensusPeaks.bed  | sort -k1,1 -k2,2n | bedtools merge > ${outdir}/merged.bed

bedtools subtract -A -a ${outdir}/merged.bed \
-b ${peakdir}/SDL_H2AUB.consensusPeaks.bed > ${outdir}/H2AUB_EMBRYOS_substract_SDL.bed

rm ${outdir}/merged.bed

# H3K4me3
bedtools subtract -A -a ${peakdir}/WS_K4ME3.consensusPeaks.bed \
-b ${peakdir}/MG_K4ME3.consensusPeaks.bed > ${outdir}/K4ME3_WS_substract_MG.bed

bedtools subtract -A -a ${peakdir}/MG_K4ME3.consensusPeaks.bed \
-b ${peakdir}/WS_K4ME3.consensusPeaks.bed > ${outdir}/K4ME3_MG_substract_WS.bed

cat ${peakdir}/WS_K4ME3.consensusPeaks.bed \
 ${peakdir}/MG_K4ME3.consensusPeaks.bed  | sort -k1,1 -k2,2n | bedtools merge > ${outdir}/merged2.bed

bedtools subtract -A -a ${outdir}/merged2.bed \
-b ${peakdir}/SDL_K4ME3.consensusPeaks.bed > ${outdir}/K4ME3_EMBRYOS_substract_SDL.bed

rm ${outdir}/merged2.bed

# H3K27me3
bedtools subtract -A -a ${peakdir}/WS_K27ME3.consensusPeaks.bed \
-b ${peakdir}/MG_K27ME3.consensusPeaks.bed > ${outdir}/K27ME3_WS_substract_MG.bed

bedtools subtract -A -a ${peakdir}/MG_K27ME3.consensusPeaks.bed \
-b ${peakdir}/WS_K27ME3.consensusPeaks.bed > ${outdir}/K27ME3_MG_substract_WS.bed


cat ${peakdir}/WS_K27ME3.consensusPeaks.bed ${peakdir}/MG_K27ME3.consensusPeaks.bed  | sort -k1,1 -k2,2n | bedtools merge > ${outdir}/merged3.bed

bedtools subtract -A -a ${outdir}/merged3.bed \
-b ${peakdir}/SDL14D_K27ME3.consensusPeaks.bed > ${outdir}/K27ME3_EMBRYOS_substract_SDL.bed

rm ${outdir}/merged3.bed