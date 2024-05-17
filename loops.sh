#!/bin/bash
#SBATCH --account=sam77_h
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=150GB
#SBATCH --partition=sla-prio
#SBATCH --job-name=MCAI54-R1-call_loops_using_mustache_normalized_cool
#SBATCH --output=../output/slurm-%x-%j.out
umask 007

source ~/conda_init.sh
conda activate mustache

SECONDS=0

~/group/lab/siddharth/software/mustache/mustache/mustache.py -f ../output/DE_merged_bio_rep_1000_normalized.cool -ch chr1 -r 1kb -st 0.7 -pt 0.05 -p 4 -o ../output/loops/loops_DE_chr1_1000_normalized.tsv

ELAPSED="Elapsed(trainNN.sh): $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo ""
echo $ELAPSED

