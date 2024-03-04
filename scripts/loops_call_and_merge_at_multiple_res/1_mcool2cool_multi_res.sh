#!/bin/bash
#SBATCH --account=sam77_h
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=150GB
#SBATCH --partition=sla-prio
#SBATCH --job-name=MCAI64-R1-mcool2cool
#SBATCH --output=../../output/slurm-%x-%j.out
umask 007

source ~/conda_init.sh
conda activate hicexplorer

SECONDS=0

in_mcool_file_path="../../input/DE_merged_bio_rep.mcool"       # input .mcool file
in_resolution=1000                                             # high resolution present in .mcool
out_resolution_list=(500000 600000 700000 800000)              # list of output resolutions (must be multiples of in_resolution)
out_path="../../output"                                        # path to output directory
out_cool_file_prefix="DE_merged_bio_rep"                       # output .cool files prefix

_in_mcool_file_path="${in_mcool_file_path}::/resolutions/${in_resolution}"


for out_resolution in "${out_resolution_list[@]}"
do
    echo "Resolution = ${out_resolution}"
unnormalized_out_cool_file_path="${out_path}/${out_cool_file_prefix}_${out_resolution}_unnormalized.cool"
    normalized_out_cool_file_path="${out_path}/${out_cool_file_prefix}_${out_resolution}_normalized.cool"

    k=$((out_resolution / in_resolution))

    echo "Factor = $k"
    echo ${_in_mcool_file_path}
    echo ${unnormalized_out_cool_file_path}
    echo ${normalized_out_cool_file_path}

    cooler coarsen -k $k -n 4 ${_in_mcool_file_path} -o ${unnormalized_out_cool_file_path}

    cp ${unnormalized_out_cool_file_path} ${normalized_out_cool_file_path}

    #cooler balance ${normalized_out_cool_file_path}

    echo "#### Normalize cool file ####"
    python normalize_cool.py ${normalized_out_cool_file_path}
done

ELAPSED="Elapsed(trainNN.sh): $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo ""
echo $ELAPSED

