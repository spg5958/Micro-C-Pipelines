#!/bin/bash
#SBATCH --account=sam77_h
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=100GB
#SBATCH --partition=sla-prio
#SBATCH --job-name=MCAI53-R1-run_pileup_analysis_PRDM1_vs_H2AK119Ub1-H2AK119Ub1.no_PRDM1
#SBATCH --output=../../output/slurm-%x-%j.out
umask 007

source ~/conda_init.sh
conda activate coolpup

SECONDS=0

echo "pileup_analysis"
python pileup_analysis.py

ELAPSED="Elapsed(trainNN.sh): $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo ""
echo $ELAPSED
