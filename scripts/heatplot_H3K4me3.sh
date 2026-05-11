#!/bin/bash
#SBATCH --job-name=heatplot
#SBATCH --output=../logs/SLURM_out/%A_%a.out
#SBATCH --error=../logs/SLURM_err/%A_%a.err
#SBATCH --array=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4G
#SBATCH --time=4:00:00

eval "$(conda shell.bash hook)"
conda activate deepTools

outdir="../results/heatplot_H3K4me3"
genesbed2000="${outdir}/TAIR10_GFF3_genes.random2000.bed"
genesbed5000="${outdir}/TAIR10_GFF3_genes.random5000.bed"
genesbed="${outdir}/TAIR10_GFF3_genes.bed"

indir="../results/heatplot_H3K4me3"
input=$(ls ${indir}/*.bw)

name="heatplots_aroundGeneBodies_plasma"
output="${outdir}/${name}"

#
computeMatrixAndPlot(){
    # $1 = input bed
    # $2 = name (2000, 5000, or all)
    computeMatrix scale-regions -S ${input} \
    -R ${1} --beforeRegionStartLength 1500 --regionBodyLength 2500 --afterRegionStartLength 1500 --skipZeros \
    --missingDataAsZero \
    -o "${output}.${2}.matrix_gene.mat.gz" -p 12 &&
    plotHeatmap -m "${output}.${2}.matrix_gene.mat.gz" -out "${output}.${2}.png" --sortUsing sum \
    --colorMap plasma \
    --heatmapHeight 20 --heatmapWidth 15
}

computeMatrixAndPlot "${genesbed2000}" "2000" &&
computeMatrixAndPlot "${genesbed5000}" "5000" &&
computeMatrixAndPlot "${genesbed}" "all"