#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate bedtools

peakdir="../results/macs2"
readsdir="../results/mergedBam_bed"
outdir="../results/peakStats"
mkdir -p "${outdir}"

echo -e "sample\tFRiP" > "${outdir}"/FRiP.tsv

for sample in $(ls ${peakdir}); do
peakfile=${peakdir}/${sample}/${sample}.broad.q0.05_peaks.broadPeak
readsfile=${readsdir}/${sample}_picard.dupRemove.bed
fragNotinPeaks=$(bedtools intersect -v -a "${readsfile}" \
-b "${peakfile}" | wc -l)
fragTotal=$(cat "${readsfile}" | wc -l)
fraginPeaks=$((fragTotal-fragNotinPeaks))
FRiP=$(echo "scale=4; ${fraginPeaks}/${fragTotal}*100" | bc)
echo -e "${sample}\t${FRiP}" >> "${outdir}"/FRiP.tsv
done