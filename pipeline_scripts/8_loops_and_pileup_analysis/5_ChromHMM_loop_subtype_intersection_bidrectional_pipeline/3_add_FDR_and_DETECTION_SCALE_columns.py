import pandas as pd
import numpy as np


# path to merged loops .bedpe file (this can be any merged loops file eg. chip_seq_and_loop_subtypes_intersection_left.bedpe)
merged_loops_df=pd.read_csv("../input/MCAI74-R1-output/chip_seq_and_loop_subtypes_intersection/chip_seq_and_loop_subtypes_intersection_left.bedpe", sep="\t", header=None)

# path to individual resolution loops.tsv file from mustache
loops_file_path = "../input/MCAI70-R1-input/mustache_loops_all_chr_multi_res/loops_chr_all_{resolution}.tsv" 

# list of resolutions
resolution_list = [1000, 2000, 4000, 8000, 10000]

# FRD filter threshold
FDR = 0.5

# path to output directory
out_path= "../output"

# output file prefix (<out_prefix>.tsv & <out_prefix>.bedpe)
out_file_prefix = "chip_seq_and_loop_subtypes_intersection_left_extended"


loops_multi_res_df_list=[]
for resolution in resolution_list:
    print(f"Resolution = {resolution}")
   
    _loops_file_path=loops_file_path.format(resolution=resolution)
    print(_loops_file_path)
    loops_df=pd.read_csv(_loops_file_path, sep="\t")
    print(loops_df)
    print(loops_df.shape)
    loops_multi_res_df_list.append(loops_df)


print("#"*20)
loops_multi_res_df_concat=pd.concat(loops_multi_res_df_list)
print(loops_multi_res_df_concat)
print(loops_multi_res_df_concat.shape)
print("#"*20)

print(merged_loops_df)
print(merged_loops_df.shape)
merged_loops_df["res1"]=merged_loops_df[2]-merged_loops_df[1]
merged_loops_df["res2"]=merged_loops_df[5]-merged_loops_df[4]
merged_loops_df["equal"]=(merged_loops_df["res1"]==merged_loops_df["res2"])

print(merged_loops_df)
print(merged_loops_df["res1"].value_counts())
print(merged_loops_df["res2"].value_counts())
print(merged_loops_df["equal"].value_counts(normalize=True))


merged_loops_df.columns=list(loops_multi_res_df_concat.columns[:6])+list(merged_loops_df.columns[6:])
print(merged_loops_df)

merged_loops_df_ext = merged_loops_df.copy()
merged_loops_df_ext=pd.merge(merged_loops_df_ext,loops_multi_res_df_concat,on=list(loops_multi_res_df_concat.columns[:6]),how="left")
print(merged_loops_df_ext)

merged_loops_df_ext_filtered=merged_loops_df_ext[merged_loops_df_ext["FDR"]<=FDR]
print(merged_loops_df_ext_filtered)
merged_loops_df_ext_filtered=merged_loops_df_ext_filtered.drop(["res1","res2","equal"],axis=1)
print(merged_loops_df_ext_filtered)

merged_loops_df_ext_filtered.to_csv(f"{out_path}/{out_file_prefix}.tsv",index=False,sep="\t")
merged_loops_df_ext_filtered.to_csv(f"{out_path}/{out_file_prefix}.bedpe", header=False, index=False, sep="\t")
