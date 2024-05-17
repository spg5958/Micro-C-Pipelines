#!/bin/bash
#SBATCH --account=sam77_h
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=150GB
#SBATCH --partition=sla-prio
#SBATCH --job-name=MCAI46-R2-3_normalize_cool
#SBATCH --output=../output/slurm-%x-%j.out
umask 007

source ~/conda_init.sh
conda activate hicexplorer

SECONDS=0

in_file_path="../output/cool_matrix/hPSC_merged_bio_rep/1000/hPSC_merged_bio_rep_1000_random_chr_removed_unnormalized.cool"
in_resolution=1000
out_path="../output/cool_matrix/hPSC_merged_bio_rep/1000"
#in_file_prefix="$(basename "${in_cool_path%.*}")"
out_file_path="${out_path}/hPSC_merged_bio_rep_1000_random_chr_removed_normalized.cool"


cp ${in_file_path} ${out_file_path}

echo "#### Normalize cool file ####"
python normalize_cool.py ${out_file_path}

ELAPSED="Elapsed(trainNN.sh): $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo ""
echo $ELAPSED

