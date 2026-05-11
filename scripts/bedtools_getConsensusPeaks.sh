#!/bin/bash
#SBATCH --job-name=bedtools_consensusPeaks
#SBATCH --output=../logs/SLURM_out/%A_%a.out
#SBATCH --error=../logs/SLURM_err/%A_%a.err
#SBATCH --array=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --time=4:00:00

# Select peaks found in > 50% of total replicates
# In the context of three replicates:
# Peaks have overlap in at least one another replicate considered reproducible peaks

eval "$(conda shell.bash hook)"
conda activate bedtools

outdir="../results/bedtools_consensusPeak"
mkdir -p ${outdir}


hist="H3K27ME3"
rep1="406231_28-JS_WT_H3k27me3_rep1_S3"
rep2="406231_34-JS_WT_H3k27me3_rep2_S27"
rep3="406231_38-JS_WT_H3k27me3_rep3_S10"

peakDir="../data_fromOtherProjects/bed_broadpeaks_new_datasets"

cat ${peakDir}/${rep1}.broad.q0.05_peaks.broadPeak \
${peakDir}/${rep2}.broad.q0.05_peaks.broadPeak \
${peakDir}/${rep3}.broad.q0.05_peaks.broadPeak \
| sort -k1,1 -k2,2n | bedtools merge -d 0 -c 1 -o count | awk '$4 > 1' > ${outdir}/"${hist}".consensusPeaks.bed

# 'bedtools merge -c 1 -o count' counts how many original peaks are included in a merged peak. This number can exceed 3
# because there could be more than one peaks in a single bed file involved in the final merged peaks
# In the case of finding peaks reproduced at least one time, there is no problem
# but for majority voting, multiIntersectBed in bedtools is needed