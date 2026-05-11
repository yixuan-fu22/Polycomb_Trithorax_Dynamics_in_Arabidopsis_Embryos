#!/bin/bash

peakDir="../results/macs2"
outdir="../results/peakStats"
mkdir -p ${outdir}
echo -e "sample\tpeakNumber\taveragePeakLength\tpercentage_genome_covered" > ${outdir}/peakNumberLength.tsv


for i in $(ls ${peakDir}); do
  file=${peakDir}/${i}/${i}.broad.q0.05_peaks.broadPeak
  line_count=$(cat "$file" | wc -l)
  avg_Length=$(cat $file | awk -F '\t' 'function abs(x) {return ((x < 0.0) ? -x : x)} {print abs($3 - $2)}' | awk '{ sum += $1; count++ } END { if (count > 0) print sum / count }')
  percent_genome_covered=$(echo "scale=6; $line_count * $avg_Length / 119481543 * 100" | bc)
  echo -e "${i}\t$line_count\t$avg_Length\t$percent_genome_covered" >> ../results/peakStats/peakNumberLength.tsv
done