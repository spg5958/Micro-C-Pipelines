#!/bin/bash
#SBATCH --account=sam77_h
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=150GB
#SBATCH --partition=sla-prio
#SBATCH --job-name=MCAI56-R1-mcool2cool
#SBATCH --output=../output/slurm-%x-%j.out
umask 007

source ~/conda_init.sh
conda activate hicexplorer

SECONDS=0

in_mcool_file_path="../input/MCAI50-R1-output/mcool_matrix/DE_merged_bio_rep/DE_merged_bio_rep.mcool"
in_resolution=1000
out_resolution=500000
out_path="../output/mcool2cool"
out_cool_file_prefix="DE_merged_bio_rep"

_in_mcool_file_path="${in_mcool_file_path}::/resolutions/${in_resolution}"
unnormalized_out_cool_file_path="${out_path}/${out_cool_file_prefix}_${out_resolution}_unnormalized.cool"
normalized_out_cool_file_path="${out_path}/${out_cool_file_prefix}_${out_resolution}_normalized.cool"

k=$((out_resolution / in_resolution))

echo "Factor = $k"
echo ${_in_mcool_file_path}
echo ${unnormalized_out_cool_file_path}
echo ${normalized_out_cool_file_path}


if [[ "$k" -eq 1 ]];
    then
        echo "k==1 using HiCExplorer"
        hicConvertFormat -m ${_in_mcool_file_path} --inputFormat cool --outputFormat cool -o ${unnormalized_out_cool_file_path}
    else
	cooler coarsen -k $k -n 4 ${_in_mcool_file_path} -o ${unnormalized_out_cool_file_path}
    fi


cp ${unnormalized_out_cool_file_path} ${normalized_out_cool_file_path}

#cooler balance ${normalized_out_cool_file_path}

echo "#### Normalize cool file ####"
python normalize_cool.py ${normalized_out_cool_file_path}

ELAPSED="Elapsed(trainNN.sh): $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo ""
echo $ELAPSED

