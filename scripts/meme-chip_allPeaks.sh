#!/bin/bash
#SBATCH --job-name=meme-chip
#SBATCH --output=../logs/SLURM_out/%A_%a.out
#SBATCH --output=../logs/SLURM_err/%A_%a.err
#SBATCH --array=0-10
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=28800
#SBATCH --time=4:00:00

# environment
eval "$(conda shell.bash hook)"
conda activate meme

# input FASTA
fasta_folder="../results/fastas_for_motif_analysis_allpeaks"
fasta=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${fasta_folder}/list.fasta.txt)

# database
database="../reference_data/meme_db/motif_databases/JASPAR/JASPAR2022_CORE_plants_non-redundant_v2.meme"

# outdir
outdir="../results/meme-chip_allPeaks_newDatasets"
mkdir -p ${outdir}

# individual results folder
x=${fasta/".consensusPeaks.fasta"/} # a REG extract string before the first . # x is an intermediate variable
clusterDir=${x}
mkdir -p ${outdir}/${clusterDir}

# ccut =	The maximum length of a sequence to use before it is trimmed to a central region of this size. A value of 0 indicates that sequences should not be trimmed.
# 
meme-chip -oc ${outdir}/${clusterDir} -time 240 -ccut 0 -dna -order 2 -minw 6 -maxw 15 -db ${database} -meme-mod anr -meme-nmotifs 10 -meme-searchsize 100000 -streme-pvt 0.05 -streme-align \
center -streme-totallength 4000000 -centrimo-score 5.0 -centrimo-ethresh 10.0 ${fasta_folder}/${fasta}
