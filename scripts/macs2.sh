#!/bin/bash
#SBATCH --job-name=macs2
#SBATCH --output=../logs/SLURM_out/%A_%a.out
#SBATCH --error=../logs/SLURM_err/%A_%a.err
#SBATCH --array=0-3
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --time=4:00:00

eval "$(conda shell.bash hook)"
conda activate MACS2


samples=(
  "H2AUB_MERGED1"
  "H2AUB_MERGED2"
  "H2AUB_MERGED3"
)

# Corresponding controls
controls=(
  "395141_07-SDL_H3_REP2_S10"
  "406231_13-YF_WT_H3_rep1_S40"
  "406231_15-YF_WT_H3_rep2_S9"
)

sample=${samples[$SLURM_ARRAY_TASK_ID]}
control=${controls[$SLURM_ARRAY_TASK_ID]}

#datatype=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${metafile}" | awk '{print $7}' | tr -d '[:space:]')
sampledir="../results/mergedBam"
controldir="../data_fromOtherProjects/bam_duplicationRemoved_new_datasets"
outdir="../results/macs2"
logdir="../logs/macs2_log"
mkdir -p "$outdir"
mkdir -p "$logdir"
q_value=0.05


(macs2 callpeak -t "$sampledir"/${sample}_picard.dupRemove.bam \
      -c "$controldir"/${control}_picard.dupRemove.bam \
      -f BAMPE -g 119481543 -n ${sample}.broad.q"$q_value" \
      --outdir "$outdir"/"$sample" \
      --nomodel \
      -q "$q_value" --broad --broad-cutoff 0.1 -B --SPMR) \
      2> "$logdir"/${sample}.broad.q"$q_value".summary.txt