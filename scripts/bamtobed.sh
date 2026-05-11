#!/bin/bash
#SBATCH --job-name=bamtobed
#SBATCH --output=../logs/SLURM_out/%A_%a.out
#SBATCH --error=../logs/SLURM_err/%A_%a.err
#SBATCH --array=0-3
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --time=4:00:00

# This scripts includes samtools_sortByName.sh so it doesn't need to be run seperately

indir="../results/mergedBam"
beddir="../results/mergedBam_bed"

samples=(
  "H2AUB_MERGED1"
  "H2AUB_MERGED2"
  "H2AUB_MERGED3"
)

sample=${samples[$SLURM_ARRAY_TASK_ID]}
mkdir -p ${beddir}

# make an intermediate dir for sorted bam
sortdir="../results/mergedBam_sortByName"
mkdir -p ${sortdir}

eval "$(conda shell.bash hook)"
conda activate samtools
samtools sort -n --output-fmt BAM ${indir}/${sample}_picard.dupRemove.bam \
-o ${sortdir}/${sample}_picard.dupRemove.sorted.bam

conda activate bedtools

bedtools bamtobed -i "${sortdir}"/"${sample}"_picard.dupRemove.sorted.bam -bedpe > "${beddir}"/"${sample}"_picard.dupRemove.bed