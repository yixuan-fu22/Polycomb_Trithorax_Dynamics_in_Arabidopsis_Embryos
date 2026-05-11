#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate bedtools

peakdir="../data_fromOtherProjects/bedtools_consensusPeak"
outdir="../results/bedtools_peaks_H3K27me3_intersect_H3K4me3"
mkdir -p ${outdir}

# Get H2Aub peaks that overlap with H3K27me3
# WS
bedtools intersect -wa -a ${peakdir}/WS_K27ME3.consensusPeaks.bed  \
-b ${peakdir}/WS_K4ME3.consensusPeaks.bed > ${outdir}/WS_H3K27ME3_OVERLAPWITH_H3K4ME3.bed

# MG
bedtools intersect -wa -a ${peakdir}/MG_K27ME3.consensusPeaks.bed \
-b ${peakdir}/MG_K4ME3.consensusPeaks.bed > ${outdir}/MG_H3K27ME3_OVERLAPWITH_H3K4ME3.bed

# SDL
bedtools intersect -wa -a ${peakdir}/SDL_H3K27ME3_newDataset.consensusPeaks.bed \
-b ${peakdir}/SDL_K4ME3.consensusPeaks.bed > ${outdir}/SDL_H3K27ME3_OVERLAPWITH_H3K4ME3.bed

