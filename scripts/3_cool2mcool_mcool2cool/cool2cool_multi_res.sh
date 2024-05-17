#!/bin/bash
#SBATCH --account=sam77_h
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=100GB
#SBATCH --partition=sla-prio
#SBATCH --job-name=MCAI47-R2-2_cool2cool_multi_res
#SBATCH --output=../output/slurm-%x.out
umask 007

source ~/conda_init.sh
conda activate hicexplorer

SECONDS=0

# input .mcool file
in_cool_file_path="../output/cool_matrix/hFG_Seq021_merged_rep/1000/hicpro2higlass_out_unnormalized/hFG_Seq021_merged_rep_1000.cool"

# .cool file resolution
in_resolution=1000

# list of output resolutions (must be multiples of in_resolution)
out_resolution_list=(1000 2000 4000 8000 10000 20000)

# path to output directory
out_path="../output/cool_files"

# output .cool files prefix
out_cool_file_prefix="hFG_Seq021_merged_rep"

mkdir -p ${out_path}

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

    if [[ "$k" -eq 1 ]]; 
    then
    	echo "k==1 copying input .cool file"
	cp ${in_cool_file_path} ${unnormalized_out_cool_file_path} 
    else	        
    	cooler coarsen -k $k -n 4 ${in_cool_file_path} -o ${unnormalized_out_cool_file_path}
    fi

    cp ${unnormalized_out_cool_file_path} ${normalized_out_cool_file_path}

    #cooler balance ${normalized_out_cool_file_path}

    echo "#### Normalize cool file ####"
    python normalize_cool.py ${normalized_out_cool_file_path}
done

ELAPSED="Elapsed(trainNN.sh): $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo ""
echo $ELAPSED

