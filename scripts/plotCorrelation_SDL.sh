#!/bin/bash
#SBATCH --job-name=plotCorrelation
#SBATCH --output=../logs/SLURM_out/%A_%a.out
#SBATCH --error=../logs/SLURM_err/%A_%a.err
#SBATCH --array=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=4:00:00

eval "$(conda shell.bash hook)"
conda activate deepTools

dir="../results/correlation_SDL_newH2AUB"

multiBigwigSummary bins --binSize 500 \
 --bwfiles $(ls ${dir}/*.bw) -p 8 -o ${dir}/multibwSummary.bins500bp.npz &&
plotCorrelation \
    -in ${dir}/multibwSummary.bins500bp.npz \
    --corMethod spearman --skipZeros \
    --plotTitle "Spearman Correlation of Read Counts" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o ${dir}/heatmap_SpearmanCorr_bins500bp.png   \
    --outFileCorMatrix ${dir}/SpearmanCorr_bins500bp.tab

# plot with out numbers # in the folder
'''
 dir="." &&   plotCorrelation \
    -in ${dir}/multibwSummary.bins500bp.npz \
    --corMethod spearman --skipZeros \
    --plotTitle "Spearman Correlation of Read Counts" \
    --whatToPlot heatmap --colorMap RdYlBu \
    -o ${dir}/heatmap_SpearmanCorr_bins500bp.noNumbers.png   \
    --outFileCorMatrix ${dir}/SpearmanCorr_bins500bp.tad
'''