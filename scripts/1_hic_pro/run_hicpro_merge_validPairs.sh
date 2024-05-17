#!/bin/bash
#SBATCH --account=sam77_h
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=100GB
#SBATCH --partition=sla-prio
#SBATCH --job-name=MCAI29-R1-run_hicpro_merged_rep
#SBATCH --output=../output/slurm-%x-%j.out
umask 007

source ~/conda_init.sh
conda activate hicpro

SECONDS=0

/storage/home/spg5958/group/lab/siddharth/software/HiC-Pro/HiC-Pro-repo/bin/HiC-Pro -i ../input/merged_validPairs_from_hic-pro_out_individual_rep-Iwafuchi_Seq021_MicroC -o ../output/hic-pro_out_merged_rep-Iwafuchi_Seq021_MicroC -s merge_persample -s build_contact_maps -s ice_norm -c config-hicpro-merged_rep.txt

ELAPSED="Elapsed(trainNN.sh): $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo ""
echo $ELAPSED
