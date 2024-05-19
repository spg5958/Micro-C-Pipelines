import numpy as np
import pandas as pd
import pybedtools
from subprocess import call


# path to ChIP-seq .bed file
chip_seq_path = "../input/dEN_SOX17.1rpm.bed"

chip_seq_tf_name = "SOX17"
 
# path to chromHMM domain .bed file
chromHMM_file_path="../input/cell1_8_dense_ID.bed"

# Left histone mark/TF
subtype_1 = "Promoters" 

# Right histone mark/TF
subtype_2 = "Enhancers"                                                    

# path to merged loops.bedpe file
loops_file_path = "../input/loops_merged_all_chr_multi_res.bedpe"   

# output path
out_path = f"../output/script_3/{chip_seq_tf_name}_{subtype_1}_{subtype_2}"     


call(["mkdir", "-p", out_path])

print("chip_seq_bed_df")
chip_seq_bed_df = pd.read_csv(chip_seq_path, sep="\t", header=None)
print(chip_seq_bed_df)
print(chip_seq_bed_df.shape)
print(chip_seq_bed_df.columns)

# extract beds
chromHMM_df = pd.read_csv(chromHMM_file_path, sep="\t")
print(chromHMM_df.shape)
print(chromHMM_df)
print(chromHMM_df.columns)
print(chromHMM_df["ChromHMM_Category"].unique())

print("loops_df")
print(loops_file_path)
loops_df=pd.read_csv(loops_file_path, sep="\t", header=None)
loops_df["loop_number"]=np.arange(loops_df.shape[0])
print(loops_df)
print(loops_df.shape)

chromHMM_bed_subtype_1_df = chromHMM_df[chromHMM_df["ChromHMM_Category"]==subtype_1]
chromHMM_bed_subtype_1_df = chromHMM_bed_subtype_1_df.iloc[:,:6]
print(chromHMM_bed_subtype_1_df)
chromHMM_bed_subtype_2_df = chromHMM_df[chromHMM_df["ChromHMM_Category"]==subtype_2]
chromHMM_bed_subtype_2_df = chromHMM_bed_subtype_2_df.iloc[:,:6]
print(chromHMM_bed_subtype_2_df)

# intersect ChIP-seq and chromHMM .bed file
chip_seq_bed = pybedtools.BedTool.from_dataframe(chip_seq_bed_df)
chromHMM_subtype_1_bed = pybedtools.BedTool.from_dataframe(chromHMM_bed_subtype_1_df)
chromHMM_subtype_2_bed = pybedtools.BedTool.from_dataframe(chromHMM_bed_subtype_2_df)
bed_1_df = chromHMM_subtype_1_bed.intersect(chip_seq_bed, u=True).to_dataframe()
bed_2_df = chromHMM_subtype_2_bed.intersect(chip_seq_bed, u=True).to_dataframe()
print("bed_1_df")
print(bed_1_df)
print(bed_1_df.shape)
print("bed_2_df")
print(bed_2_df)
print(bed_2_df.shape)

def intersect(loops_df, bed_1_df, bed_2_df):
    
    # intersection
    loops_df.columns = ["BIN1_CHR","BIN1_START","BIN1_END","BIN2_CHROMOSOME","BIN2_START","BIN2_END","loop_number"]
    
    loops_left_df=loops_df[["BIN1_CHR","BIN1_START","BIN1_END","loop_number"]]

    loops_right_df=loops_df[["BIN2_CHROMOSOME","BIN2_START","BIN2_END","loop_number"]]

    loops_left_bed=pybedtools.BedTool.from_dataframe(loops_left_df)

    loops_right_bed=pybedtools.BedTool.from_dataframe(loops_right_df)

    # https://github.com/daler/pybedtools/issues/238
    bed_1 = pybedtools.BedTool.from_dataframe(bed_1_df)
    bed_2 = pybedtools.BedTool.from_dataframe(bed_2_df)

    loops_left_bed_1_intersect=loops_left_bed.intersect(bed_1, u=True).to_dataframe()
    loops_left_bed_1_intersect=loops_left_bed_1_intersect.rename(columns={"chrom":"BIN1_CHR", "start":"BIN1_START", "end":"BIN1_END", "name":"loop_number"})

    loops_right_bed_2_intersect=loops_right_bed.intersect(bed_2, u=True).to_dataframe()
    loops_right_bed_2_intersect=loops_right_bed_2_intersect.rename(columns={"chrom":"BIN2_CHROMOSOME", "start":"BIN2_START", "end":"BIN2_END", "name":"loop_number"})

    loops_bed_1_on_left=loops_df.merge(loops_left_bed_1_intersect, how="inner", on=["BIN1_CHR", "BIN1_START", "BIN1_END", "loop_number"])

    loops_bed_2_on_right = loops_df.merge(loops_right_bed_2_intersect, how="inner", on=["BIN2_CHROMOSOME", "BIN2_START", "BIN2_END", "loop_number"])

    loops_bed_1_on_left_bed_2_on_right = loops_bed_1_on_left.merge(loops_bed_2_on_right, how="inner")
    loops_bed_1_on_left_bed_2_on_right = loops_bed_1_on_left_bed_2_on_right.iloc[:,:6]
    print(loops_bed_1_on_left_bed_2_on_right)
    print(loops_bed_1_on_left_bed_2_on_right.shape)
    
    # verify
    verify_df = loops_bed_2_on_right.merge(loops_bed_1_on_left, how="inner")
    if verify_df.shape[0] == loops_bed_1_on_left_bed_2_on_right.shape[0]:
        print(f"{verify_df.shape[0]} == {loops_bed_1_on_left_bed_2_on_right.shape[0]}")
    else:
        print(f"Error: {verify_df.shape[0]} != {loops_bed_1_on_left_bed_2_on_right.shape[0]}")
        raise Exception("an error occurred")
        
    return loops_bed_1_on_left_bed_2_on_right


for left_mark, right_mark in [(subtype_1, subtype_1), (subtype_2, subtype_2), (subtype_1, subtype_2)]:
    
    if left_mark!=right_mark:
        print("Left and Right subtypes are not equal")
        
         # path to output file
        out_file_path_1 = f"{out_path}/loops_{left_mark}_{chip_seq_tf_name}_left_{right_mark}_{chip_seq_tf_name}_right.bedpe"
        # path to output file
        out_file_path_2 = f"{out_path}/loops_{right_mark}_{chip_seq_tf_name}_left_{left_mark}_{chip_seq_tf_name}_right.bedpe"        
        
        loops_bed_1_on_left_bed_2_on_right = intersect(loops_df, bed_1_df, bed_2_df)
        loops_bed_2_on_left_bed_1_on_right = intersect(loops_df, bed_2_df, bed_1_df)

        loops_bed_1_on_left_bed_2_on_right.to_csv(out_file_path_1, sep="\t", header=None, index=False)
        loops_bed_2_on_left_bed_1_on_right.to_csv(out_file_path_2, sep="\t", header=None, index=False)
    elif (left_mark==right_mark) and (left_mark==subtype_1):
        print(f"Left and Right subtypes are equal - {left_mark}")

        # path to output file
        out_file_path = f"{out_path}/loops_{left_mark}_{chip_seq_tf_name}_left_{right_mark}_{chip_seq_tf_name}_right.bedpe"        
        
        loops_bed_1_on_left_bed_2_on_right = intersect(loops_df, bed_1_df, bed_1_df)
     
        loops_bed_1_on_left_bed_2_on_right.to_csv(out_file_path, sep="\t", header=None, index=False)
    elif (left_mark==right_mark) and (left_mark==subtype_2):
        print(f"Left and Right subtypes are equal - {left_mark}")

        # path to output file
        out_file_path = f"{out_path}/loops_{left_mark}_{chip_seq_tf_name}_left_{right_mark}_{chip_seq_tf_name}_right.bedpe"        
        
        loops_bed_1_on_left_bed_2_on_right = intersect(loops_df, bed_2_df, bed_2_df)
     
        loops_bed_1_on_left_bed_2_on_right.to_csv(out_file_path, sep="\t", header=None, index=False)
