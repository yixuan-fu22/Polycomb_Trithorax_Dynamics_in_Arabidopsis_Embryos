#!/bin/bash
#SBATCH --job-name=samtoolsMergeRep
#SBATCH --output=../logs/SLURM_out/%A_%a.out
#SBATCH --error=../logs/SLURM_err/%A_%a.err
#SBATCH --array=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=4:00:00

eval "$(conda shell.bash hook)"
conda activate samtools

## use sorted, filterred bam in samtools_sort to merge
indir="../data_fromOtherProjects/bam_duplicationRemoved_new_datasets"
outdir="../results/mergedBam"
mkdir -p ${outdir}

##
rep1="406231_14-YF_WT_H2Aub_rep1_S7_picard.dupRemove.bam"
rep2="406231_16-YF_WT_H2Aub_rep2_S5_picard.dupRemove.bam"
fileName="H2AUB_MERGED2_picard.dupRemove.bam"
##
samtools merge --threads 8 -f -o ${outdir}/${fileName} ${indir}/${rep1} \
${indir}/${rep2}

## count the number of reads to confirm
nreadRep1=$(samtools view --threads 8 -c ${indir}/${rep1})
nreadRep2=$(samtools view --threads 8 -c ${indir}/${rep2})
nreadMerge=$(samtools view --threads 8 -c ${outdir}/${fileName})

echo rep1 ${nreadRep1}
echo rep2 ${nreadRep2}
echo fileName ${nreadMerge}

