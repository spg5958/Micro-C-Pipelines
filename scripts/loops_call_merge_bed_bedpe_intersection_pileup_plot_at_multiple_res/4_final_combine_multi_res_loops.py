import numpy as np
import pandas as pd
import subprocess
import os
import io
import pybedtools


loops_dir="../output/merge_loops"    # path to directory containing loops at different resolutions in bedpe format
out_path="../output/merge_loops"     # path to output directory

loops_df_list=[]
for file_name in os.listdir(loops_dir):
    if file_name.endswith(".bedpe"):
        print(file_name)
        df=pd.read_csv(f"{loops_dir}/{file_name}",sep="\t",header=None)
        df=df.dropna(axis=1)
        df["file_name"]=file_name
        print(df)
        print()
        loops_df_list.append(df)
        #loops_file_path_list.append(f"{loops_dir}/{file_name}")
print(len(loops_df_list))

loops_df_orig=pd.concat(loops_df_list,ignore_index=True)
print(loops_df_orig)
print(loops_df_orig.shape)

loops_df=loops_df_orig.copy()
loops_df["status"]="nan"
for i in loops_df.index[:-1]:
    print("#"*20)
    print(f"i = {i}")
    row_a=loops_df.loc[[i],loops_df.columns[:-2]]
    a=pybedtools.BedTool.from_dataframe(row_a)
    print(a)
    compare_df=loops_df.loc[i+1:,:]#.reset_index(drop=True)
    print(len(compare_df))
    print(compare_df)
    #print()
    for j in compare_df.index:
        print(f"j = {j}")
        row_b=compare_df.loc[[j],:]
        if row_b.loc[j]["status"]=="drop":
            print(f"Skipping row = {j}")
            print(row_b)
            continue
        row_b=row_b.drop(["file_name","status"],axis=1)
        b=pybedtools.BedTool.from_dataframe(row_b)
        print(a)
        print(b)
        intersection=pybedtools.bedtool.BedTool.pair_to_pair(a,b, type="both")
        intersection=intersection.to_dataframe(header=None).dropna(axis=1).drop_duplicates()
        print(intersection)
        print(type(intersection))
        print(len(intersection))
        if len(intersection)==0:
            print("No intersection")
        else:
            print("intersection")
            print(loops_df.shape)
            #loops_df=loops_df.drop(j)
            print(loops_df.shape)
            print(row_a)
            print(row_b)
            res_a=row_a.loc[i,2]-row_a.loc[i,1]
            res_b=row_b.loc[j,2]-row_b.loc[j,1]
            print(res_a,res_b)
            if res_a<=res_b:
                loops_df.loc[i,"status"]="select"
                loops_df.loc[j,"status"]="drop"
            else:
                loops_df.loc[i,"status"]="drop"
                loops_df.loc[j,"status"]="select"
                break
        print("-"*20)
print(loops_df)
loops_df.to_csv(f"{out_path}/loops_merged.tsv", index=False, sep="\t")
bool_mask=loops_df["status"]=="drop"
loops_df=loops_df[~bool_mask]
loops_df.iloc[:,:6].to_csv(f"{out_path}/loops_merged.bedpe", header=False, index=False, sep="\t")
