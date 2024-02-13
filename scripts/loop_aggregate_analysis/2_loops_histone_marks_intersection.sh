#!/bin/bash
#SBATCH --account=sam77_h
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=50GB
#SBATCH --partition=sla-prio
#SBATCH --job-name=MCAI63-R1-2_loops_histone_marks_intersection
#SBATCH --output=../output/slurm-%x-%j.out
umask 007

source ~/conda_init.sh
conda activate coolpup

SECONDS=0

echo "loops_histone_marks_intersection"
python 2_loops_histone_marks_intersection.py

ELAPSED="Elapsed(trainNN.sh): $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo ""
echo $ELAPSED
