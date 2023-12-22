#!/bin/bash
#SBATCH --account=sam77_h
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=40GB
#SBATCH --partition=sla-prio
#SBATCH --job-name=MCAI28-R1-create_matrix_ab_comp
#SBATCH --output=../../output/slurm-%x-%j.out
umask 007

source ~/conda_init.sh
conda activate hicexplorer

RESOLUTION=500000
cool_file_path="../input/MCAI41-R1-output/hFG-lib548+lib549_merged/${RESOLUTION}/hFG-lib548+lib549_merged.cool"
out_path="../output/matrix_cool_ab_comp"

CHR="chr1"

mkdir -p ${out_path}

hicCorrectMatrix diagnostic_plot -m ${cool_file_path} -o ${out_path}/diagnostic_plot_${CHR}_${RESOLUTION}.png --chromosomes ${CHR}

hicCorrectMatrix correct --matrix ${cool_file_path} --correctionMethod ICE --chromosomes $CHR --iterNum 5000 --outFileName ${out_path}/corrected_ICE_${CHR}_${RESOLUTION}.cool --filterThreshold -1.8 2.6

hicTransform -m ${out_path}/corrected_ICE_${CHR}_${RESOLUTION}.cool --method obs_exp -o ${out_path}/obs_exp_${SAMPLE}_${CHR}_${RESOLUTION}.cool --chromosomes ${CHR}

hicTransform -m ${out_path}/obs_exp_${SAMPLE}_${CHR}_${RESOLUTION}.cool --method pearson -o ${out_path}/pearson_obs_exp_${SAMPLE}_${CHR}_${RESOLUTION}.cool --chromosomes ${CHR}

hicPCA -m ${out_path}/pearson_obs_exp_${SAMPLE}_${CHR}_${RESOLUTION}.cool -o ${out_path}/pca1_${SAMPLE}_${CHR}_${RESOLUTION}.bedgraph --format bedgraph -we 1 --chromosomes ${CHR}

