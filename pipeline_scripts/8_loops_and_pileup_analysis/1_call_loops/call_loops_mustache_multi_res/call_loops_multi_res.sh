#!/bin/bash
#SBATCH --account=open
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=50GB
#SBATCH --partition=open
#SBATCH --job-name=MCAI67-R1-2_call_loops_multi_res
#SBATCH --output=../output/slurm-%x-%j.out
umask 007

source ~/conda_init.sh
conda activate mustache

SECONDS=0

python call_loops_multi_res.py

ELAPSED="Elapsed(trainNN.sh): $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo ""
echo $ELAPSED

