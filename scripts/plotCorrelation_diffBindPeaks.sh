#!/bin/bash
#SBATCH --job-name=correlation
#SBATCH --output=../logs/SLURM_out/%A_%a.out
#SBATCH --output=../logs/SLURM_err/%A_%a.err
#SBATCH --array=0-6
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=28800
#SBATCH --time=4:00:00

eval "$(conda shell.bash hook)"
conda activate deepTools

DIRS=(
  "cr_MGvsSDL_H2AUB"
  "cr_MGvsSDL_K4ME3"
  "cr_WSvsMG_H3K27ME3"
  "cr_MGvsSDL_H3K27ME3"
  "cr_WSvsMG_H2AUB"
  "cr_WSvsMG_K4ME3"
)

dir="../results/${DIRS[$SLURM_ARRAY_TASK_ID]}"

id="diffBindPeaks"


 # dir="../results/cr_MGvsSDL_H2AUB"
 # dir="../results/cr_MGvsSDL_K4ME3"
 # dir="../results/cr_WSvsMG_H3K27ME3"
 # dir="../results/cr_MGvsSDL_H3K27ME3"
 # dir="../results/cr_WSvsMG_H2AUB"
  dir="../results/cr_WSvsMG_K4ME3"

multiBigwigSummary BED-file --BED $(ls ${dir}/*.bed) \
 --bwfiles $(ls ${dir}/*.bw) -p 8 -o ${dir}/multibwSummary.${id}.npz &&
plotCorrelation \
    -in ${dir}/multibwSummary.${id}.npz \
    --corMethod spearman --skipZeros \
    --plotTitle "Spearman Correlation of Read Counts" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o ${dir}/heatmap_SpearmanCorr_${id}.png   \
    --outFileCorMatrix ${dir}/SpearmanCorr_${id}.tab &&
plotPCA -in ${dir}/multibwSummary.${id}.npz \
-o ${dir}/PCA_${id}.png \
-T "PCA of read counts"

# plot with out numbers # in the folder
'''

'''

