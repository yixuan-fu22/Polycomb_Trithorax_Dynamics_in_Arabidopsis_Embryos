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


# in the folder
multiBigwigSummary bins --binSize 500 \
 --bwfiles $(ls *.bw) -p 8 -o multibwSummary.bins500bp.npz &&
plotCorrelation \
    -in multibwSummary.bins500bp.npz \
    --corMethod spearman --skipZeros \
    --plotTitle "Spearman Correlation of Read Counts" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o heatmap_SpearmanCorr_bins500bp.png   \
    --outFileCorMatrix SpearmanCorr_bins500bp.tab &&
plotPCA -in multibwSummary.bins500bp.npz \
-o PCA_bins500bp.png \
-T "PCA of read counts"

#multiBigwigSummary bins --binSize 500 \
# --bwfiles $(ls *WS*.bw *MG*.bw) -o multibwSummary.em.bins500bp.npz &&
#plotCorrelation \
#    -in multibwSummary.em.bins500bp.npz \
#    --corMethod spearman --skipZeros \
#    --plotTitle "Spearman Correlation of Read Counts" \
#    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
#    -o heatmap_SpearmanCorr_bins500bp.em.png   \
#    --outFileCorMatrix SpearmanCorr_bins500bp.em.tab &&
#plotPCA -in multibwSummary.em.bins500bp.npz \
#-o PCA_bins500bp.em.png \
#-T "PCA of read counts"