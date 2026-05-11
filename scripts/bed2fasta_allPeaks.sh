#!/bin/bash
#SBATCH --job-name=bed2fasta
#SBATCH --output=../logs/SLURM_out/%A_%a.out
#SBATCH --output=../logs/SLURM_err/%A_%a.err
#SBATCH --array=1-10
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=14400
#SBATCH --time=1:00:00

# environment
eval "$(conda shell.bash hook)"
conda activate meme

# reference
genome_file="../reference_data/refGenome/TAIR10_chr_all.fas"

# list of bed files made by (ls */*.bed > list.cluster.allpeaks.bed.txt)

# input bed
bed_folder="../data_fromOtherProjects/bedtools_consensusPeak"
bed=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${bed_folder}/list.bed.txt)

# bed2fasta
# output fasta into the same folder
echo $(bed2fasta -o ${bed_folder}/${bed/".bed"/".fasta"} ${bed_folder}/${bed} ${genome_file})

'
genome_file="../../reference_data/refGenome/TAIR10_chr_all.fas"
bed_folder="."
fasta_folder="../fastas_for_motif_analysis_allpeaks"
for bed in $(ls *.bed); do
bed2fasta -o ${fasta_folder}/${bed/".bed"/".fasta"} ${bed_folder}/${bed} ${genome_file}
done
'

