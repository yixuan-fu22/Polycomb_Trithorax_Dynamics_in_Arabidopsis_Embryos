#!/bin/bash
#SBATCH --job-name=bdg_to_bigwig
#SBATCH --output=../logs/SLURM_out/%A_%a.out
#SBATCH --error=../logs/SLURM_err/%A_%a.err
#SBATCH --array=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --time=4:00:00

eval "$(conda shell.bash hook)"
conda activate MACS2

sample="MERGED_395141_02_06_10_SDL_H2AUB_ALLREPS"
q_value=0.05
method=FE

bdgdir="../results/macs2"
bdgcmpOutdir="../results/bdgcmp"
mkdir -p ${bdgcmpOutdir}


macs2 bdgcmp -t "$bdgdir"/${sample}/${sample}.broad.q"$q_value"_treat_pileup.bdg \
-c "$bdgdir"/${sample}/${sample}.broad.q"$q_value"_control_lambda.bdg \
--outdir ${bdgcmpOutdir} --o-prefix ${sample} -m ${method}

# sortbdg

bdgsortDir="../results/bdgsort"
mkdir -p "${bdgsortDir}"

sort -k1,1 -k2,2n ${bdgcmpOutdir}/${sample}_${method}.bdg \
> "$bdgsortDir"/${sample}_${method}.sorted.bdg

# bdg to bigwig 

chromaSizeDir="../meta"
bwdir="../results/bigwig_macs2"
mkdir -p "${bwdir}"

conda activate bedGraphToBigWig

bedGraphToBigWig "${bdgsortDir}"/${sample}_${method}.sorted.bdg \
"${chromaSizeDir}"/TAIR10.chrom.size "${bwdir}"/${sample}_${method}.bw