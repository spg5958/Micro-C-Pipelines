#!/bin/bash
#SBATCH --account=sam77_h
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=50GB
#SBATCH --partition=sla-prio
#SBATCH --job-name=MCAI47-R2-1_create_cool_from_hic_pro_matrix
#SBATCH --output=../output/slurm-%x.out
umask 007

source ~/conda_init.sh
conda activate hicexplorer

in_res=1000
matrix_in_path="../input/hic_pro_out_merged_rep_hFG_Seq021/hic_results/matrix/hFG_Seq021_merged_rep/raw/1000/hFG_Seq021_merged_rep_1000.matrix"
info_file_path="../input/hg38.info"
out_path="../output/cool_matrix"

SAMPLE="hFG_Seq021_merged_rep"
echo "##############################"
echo "## ${SAMPLE} ${in_res}"
echo "##############################"

_out_path=${out_path}/${SAMPLE}/${in_res}/hicpro2higlass_out_unnormalized
mkdir -p ${_out_path}

if [ ! -f ${out_path}/${SAMPLE}/${RESOLUTION}/${SAMPLE}_${in_res}.mcool ]; then
	ls ${matrix_in_path}
	echo ${out_path}
	echo ".cool file does not exists"s
	/storage/home/spg5958/group/lab/siddharth/software/HiC-Pro/HiC-Pro-repo/bin/utils/hicpro2higlass.sh -i ${matrix_in_path} -r ${in_res} -c ${info_file_path} -o ${_out_path} -p 4 -t ${_out_path}
else
	echo "Cool file exists"
fi
