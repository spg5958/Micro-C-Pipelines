#!/bin/bash
#SBATCH --account=sam77_h
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=50GB
#SBATCH --partition=sla-prio
#SBATCH --job-name=MCAI53-R1-cool2cool
#SBATCH --output=../output/slurm-%x-%j.out
umask 007

source ~/conda_init.sh
conda activate hicexplorer

SECONDS=0

in_mcool_file_path="../input/MCAI50-R1-output/mcool_matrix/DE_merged_bio_rep/DE_merged_bio_rep.mcool::/resolutions/1000"
out_cool_file_path="../output/DE_merged_bio_rep_1000.cool"

hicConvertFormat --matrices $in_mcool_file_path --outFileName $out_cool_file_path --inputFormat cool --outputFormat cool 

ELAPSED="Elapsed(trainNN.sh): $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo ""
echo $ELAPSED

