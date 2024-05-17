#!/bin/bash
#SBATCH --account=sam77_h
#SBATCH --time=300:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=100GB
#SBATCH --partition=sla-prio
#SBATCH --job-name=MCAI29-R2-run_hicpro_merge_rep
#SBATCH --output=../output/slurm-%x-%j.out
umask 007

# activate hic-pro env
source ~/conda_init.sh
conda activate hicpro

/storage/home/spg5958/group/lab/siddharth/software/HiC-Pro/HiC-Pro-repo/bin/HiC-Pro -i ../input/merged_fastq_Seq021_MicroC -o ../output/hic_pro_out_merged_rep_hFG_Seq021 -c config-hicpro.txt
