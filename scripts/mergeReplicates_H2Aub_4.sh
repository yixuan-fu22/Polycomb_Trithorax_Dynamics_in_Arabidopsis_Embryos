#!/bin/bash
#SBATCH --job-name=samtoolsMergeRep
#SBATCH --output=../logs/SLURM_out/%A_%a.out
#SBATCH --error=../logs/SLURM_err/%A_%a.err
#SBATCH --array=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=4:00:00

# to merge H2Aub in a different way: old rep1 

eval "$(conda shell.bash hook)"
conda activate samtools

## use sorted, filterred bam in samtools_sort to merge
indir="../data_fromOtherProjects/bam_duplicationRemoved_H3K4me3_and_H2Aub/"
outdir="../results/mergedBam"
mkdir -p ${outdir}

##
rep1="395141_02-SDL_H2AUB_REP1_S22_picard.dupRemove.bam"
rep2="395141_06-SDL_H2AUB_REP2_S7_picard.dupRemove.bam"
rep3="395141_10-SDL_H2AUB_REP3_S13_picard.dupRemove.bam"
fileName="MERGED_395141_02_06_10_SDL_H2AUB_ALLREPS_picard.dupRemove.bam"
##
samtools merge --threads 8 -f -o ${outdir}/${fileName} ${indir}/${rep1} \
${indir}/${rep2} ${indir}/${rep3}

## count the number of reads to confirm
nreadRep1=$(samtools view --threads 8 -c ${indir}/${rep1})
nreadRep2=$(samtools view --threads 8 -c ${indir}/${rep2})
nreadRep3=$(samtools view --threads 8 -c ${indir}/${rep3})
nreadMerge=$(samtools view --threads 8 -c ${outdir}/${fileName})

echo rep1 ${nreadRep1}
echo rep2 ${nreadRep2}
echo rep3 ${nreadRep3}
echo fileName ${nreadMerge}

