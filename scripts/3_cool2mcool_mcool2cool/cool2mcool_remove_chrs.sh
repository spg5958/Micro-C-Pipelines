#!/bin/bash
#SBATCH --account=sam77_h
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=100GB
#SBATCH --partition=sla-prio
#SBATCH --job-name=MCAI46-R2-2_create_mcool_from_cool
#SBATCH --output=../output/slurm-%x.out
umask 007

source ~/conda_init.sh
conda activate hicexplorer

in_cool_path="../output/cool_matrix/hPSC_merged_bio_rep/1000/hicpro2higlass_out_unnormalized/hPSC_merged_bio_rep_1000.cool"
out_path="../output/cool_matrix/hPSC_merged_bio_rep/1000"

mkdir -p ${out_path}

in_cool_file_prefix="$(basename "${in_cool_path%.*}")"
out_file_prefix=${in_cool_file_prefix}_random_chr_removed_unnormalized      

echo "Prefix = ${in_cool_file_name}"

if [ ! -f ${out_path}/${out_file_prefix}.mcool ]; then
	echo ${in_cool_path}
	echo ${_out_path}
	echo "mcool file does not exists"

	echo "Remove chrs"
	hicAdjustMatrix -m ${in_cool_path} --action keep --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY -o ${out_path}/${out_file_prefix}.cool

	echo "Creating mcool"
	hicConvertFormat -m ${out_path}/${out_file_prefix}.cool --inputFormat cool --outputFormat mcool -o ${out_path}/${out_file_prefix}.mcool --resolutions 1000 2000 4000 8000 10000 20000
         
else
	echo "mcool file exists"
fi
