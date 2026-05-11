#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate bedtools

peakdir="../data_fromOtherProjects/bedtools_consensusPeak"
outdir="../results/bedtools_peaks_H2Aub_intersect_H3K27me3"
mkdir -p ${outdir}

# Get H2Aub peaks that overlap with H3K27me3
# WS
bedtools intersect -wa -a ${peakdir}/WS_H2AUB.consensusPeaks.bed \
-b ${peakdir}/WS_K27ME3.consensusPeaks.bed > ${outdir}/WS_H2AUB_OVERLAPWITH_H3K27ME3.bed

bedtools intersect -wa -v -a ${peakdir}/WS_H2AUB.consensusPeaks.bed \
-b ${peakdir}/WS_K27ME3.consensusPeaks.bed > ${outdir}/WS_H2AUB_WITHOUT_H3K27ME3.bed

# MG
bedtools intersect -wa -a ${peakdir}/MG_H2AUB.consensusPeaks.bed \
-b ${peakdir}/MG_K27ME3.consensusPeaks.bed > ${outdir}/MG_H2AUB_OVERLAPWITH_H3K27ME3.bed

bedtools intersect -wa -v -a ${peakdir}/MG_H2AUB.consensusPeaks.bed \
-b ${peakdir}/MG_K27ME3.consensusPeaks.bed > ${outdir}/MG_H2AUB_WITHOUT_H3K27ME3.bed

# SDL
bedtools intersect -wa -a ${peakdir}/MERGED_395141_02_06_10_SDL_H2AUB_ALLREPS.broad.q0.05_peaks.broadPeak.bed \
-b ${peakdir}/SDL_H3K27ME3_newDataset.consensusPeaks.bed > ${outdir}/SDL_H2AUB_OVERLAPWITH_H3K27ME3.bed

bedtools intersect -wa -v -a ${peakdir}/MERGED_395141_02_06_10_SDL_H2AUB_ALLREPS.broad.q0.05_peaks.broadPeak.bed \
-b ${peakdir}/SDL_H3K27ME3_newDataset.consensusPeaks.bed > ${outdir}/SDL_H2AUB_WITHOUT_H3K27ME3.bed


