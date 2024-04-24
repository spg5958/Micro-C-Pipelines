import numpy as np
import pandas as pd
import pybedtools
import os
import glob
from subprocess import call


# list of loop resolutions
resolution_list=[1000, 2000, 4000, 8000, 10000]      

# path to directory containing loops at different resolutions in bedpe format
loops_file_path="../output/1_loops_histone_marks_intersection_chromHMM_bidirectional/Bivalent_Bivalent/loops_left_Bivalent_right_Bivalent_{resolution}.bedpe"    

# path to output directory
out_path="../output/2_combine_multi_res_loops"    

# out file prefix
out_file_prefix="loops_merged_mutli_res_Bivalent_Bivalent"


call(["mkdir", "-p", out_path])

loops_df_list=[]
for resolution in resolution_list:
    _loops_file_path=loops_file_path.format(resolution=resolution)
    file_name=_loops_file_path.split("/")[-1]
    print(file_name)
    df=pd.read_csv(_loops_file_path,sep="\t",header=None)
    df=df.dropna(axis=1)
    df["file_name"]=file_name
    print(df)
    print()
    loops_df_list.append((file_name,df))

def sort_loops(row):
    _row=row.copy()
    if _row["end1"]<_row["start1"]:
        start1=_row["start1"].copy()
        _row["start1"]=_row["end1"]
        _row["end1"]=start1
    if _row["end2"]<_row["start2"]:
        start2=_row["start2"].copy()
        _row["start2"]=_row["end2"]
        _row["end2"]=start2
    if _row["start2"]<_row["start1"]:
        tmp=_row[["chrom1","start1","end1"]].copy()
        _row[["chrom1","start1","end1"]]=row[["chrom2","start2","end2"]]
        _row[["chrom2","start2","end2"]]=tmp
    return _row
 
prev_loops_file_name=loops_df_list[0][0]
prev_loops_df=loops_df_list[0][1]
prev_loops_df=prev_loops_df.iloc[:,:-1]
final_df=None
for i,(curr_loops_file_name, curr_loops_df) in enumerate(loops_df_list[1:]):
    print(f"prev_loops_file_name = {prev_loops_file_name}")
    print(f"curr_loops_file_name = {curr_loops_file_name}")
    a=pybedtools.BedTool.from_dataframe(prev_loops_df.iloc[:,:])
    b=pybedtools.BedTool.from_dataframe(curr_loops_df.iloc[:,:-1])
    
    a_unique=pybedtools.bedtool.BedTool.pair_to_pair(a,b, type="neither").to_dataframe(header=None).dropna(axis=1)
    if not a_unique.empty:
        a_unique.columns=["chrom1","start1","end1","chrom2","start2","end2"]#,"file_name"]
    b_unique=pybedtools.bedtool.BedTool.pair_to_pair(b,a, type="neither").to_dataframe(header=None).dropna(axis=1)
    if not b_unique.empty:
        b_unique.columns=["chrom1","start1","end1","chrom2","start2","end2"]#,"file_name"] 
    a_b_either=pybedtools.bedtool.BedTool.pair_to_pair(a,b, type="either").to_dataframe(header=None).dropna(axis=1)
    if not a_b_either.empty:
        a_either=a_b_either.iloc[:,:6]
        a_either.columns=["chrom1","start1","end1","chrom2","start2","end2"]#,"file_name"]
        b_either=a_b_either.iloc[:,6:]
        b_either.columns=["chrom1","start1","end1","chrom2","start2","end2"]#,"file_name"]
    a_b_both=pybedtools.bedtool.BedTool.pair_to_pair(a,b, type="both").to_dataframe(header=None).dropna(axis=1)
    if not a_b_both.empty:
        a_both=a_b_both.iloc[:,:6]
        a_both.columns=["chrom1","start1","end1","chrom2","start2","end2"]#,"file_name"]
        b_both=a_b_both.iloc[:,6:]
        b_both.columns=["chrom1","start1","end1","chrom2","start2","end2"]#,"file_name"]
    df_concat=pd.concat([a_unique,b_unique,a_either,b_either,a_both], axis=0)
    df_concat=df_concat.apply(sort_loops, axis=1)
    df_concat=df_concat.drop_duplicates(subset=["chrom1","start1","end1","chrom2","start2","end2"],keep='last')
    b_both=b_both.apply(sort_loops, axis=1)
    b_both=b_both.drop_duplicates(subset=["chrom1","start1","end1","chrom2","start2","end2"],keep="last")
    
    final_df=pd.merge(df_concat,b_both, indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1)
            
    verify_df=pd.concat([df_concat,b_both], axis=0)
    verify_df=verify_df.apply(sort_loops, axis=1)
    verify_df=verify_df.drop_duplicates(subset=["chrom1","start1","end1","chrom2","start2","end2"],keep=False)
    
    if final_df.shape[0]==verify_df.shape[0]:
        print(f"Verified: {final_df.shape[0]} == {verify_df.shape[0]}")
    else:
        print(f"Error during verification: {final_df.shape[0]} != {verify_df.shape[0]}")
        raise Exception("Error")
        
    prev_loops_file_name=curr_loops_file_name
    prev_loops_df=final_df.copy()

#final_df=final_df.sort_values(by=["chrom1","start1"])
print(final_df)
print(final_df.shape)
final_df.to_csv(f"{out_path}/{out_file_prefix}.tsv",index=False,sep="\t")
final_df.iloc[:,:6].to_csv(f"{out_path}/{out_file_prefix}.bedpe", header=False, index=False, sep="\t")
