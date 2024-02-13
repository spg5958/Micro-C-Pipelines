#!/bin/bash
#SBATCH --account=sam77_h
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=150GB
#SBATCH --partition=sla-prio
#SBATCH --job-name=MCAI63-R1-1_loop_calling_using_mustache
#SBATCH --output=../output/slurm-%x-%j.out
umask 007

source ~/conda_init.sh
conda activate mustache

SECONDS=0

~/group/lab/siddharth/software/mustache/mustache/mustache.py -f ../input/shared_data_Iwafuchi_lab/DE_merged_bio_rep_1000_normalized.cool -ch chr1 -r 1kb -st 0.7 -pt 0.1 -p 4 -o ../output/loops_check/loops.tsv

ELAPSED="Elapsed(trainNN.sh): $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo ""
echo $ELAPSED

