#!/bin/bash
#SBATCH --job-name=bamcoverage
#SBATCH --output=../logs/SLURM_out/%A_%a.out
#SBATCH --error=../logs/SLURM_err/%A_%a.err
#SBATCH --array=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --time=4:00:00


eval "$(conda shell.bash hook)"


indir="../results/mergedBam/"
sample="H2AUB_MERGED3_picard.dupRemove.bam"
outdir="../results/bigwig_bamcoverage"
mkdir -p ${outdir}

# first index the bam file

conda activate samtools
samtools index ${indir}/${sample} -o ${indir}/${sample}.bai

# bamcoverage

conda activate deepTools
bamCoverage -b ${indir}/${sample} -o ${outdir}/SeqDepthNorm_${sample/.bam/}.bw --binSize 10 \
--normalizeUsing RPGC --effectiveGenomeSize 119481543 --ignoreForNormalization chrC,chrM --extendReads --numberOfProcessors 4