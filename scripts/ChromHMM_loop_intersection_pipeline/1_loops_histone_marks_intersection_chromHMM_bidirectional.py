import numpy as np
import pandas as pd
import pybedtools
from subprocess import call

chromHMM_file_path = "../input/bed_files/cell1_8_dense_ID.bed"

# Left histone mark/TF
left_mark = "Polycomb" 

# Right histone mark/TF
right_mark = "Polycomb"                                                    
  
# list of resolutions
resolution_list = [1000, 2000, 4000, 8000, 10000] 

# path to loops.tsv file
loops_file_path = "../output/mustache_loops_all_chr_multi_res/loops_chr_all_{resolution}.tsv"   

# output path
out_path = "../output/1_loops_histone_marks_intersection_chromHMM_bidirectional"     


call(["mkdir", "-p", out_path])

# extract beds
chromHMM_df = pd.read_csv(chromHMM_file_path, sep="\t")
print(chromHMM_df.shape)
print(chromHMM_df)
print(chromHMM_df.columns)
print(chromHMM_df["ChromHMM_Category"])
print(chromHMM_df["ChromHMM_Category"].unique())


bed_1_df = chromHMM_df[chromHMM_df["ChromHMM_Category"]==left_mark]
bed_1_df = bed_1_df.iloc[:,:6]
print(bed_1_df)
bed_2_df = chromHMM_df[chromHMM_df["ChromHMM_Category"]==right_mark]
bed_2_df = bed_2_df.iloc[:,:6]
print(bed_2_df)


def intersect(loops_df, bed_1_df, bed_2_df):
    
    # intersection
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


for resolution in resolution_list:
    print(f"Resolution = {resolution}")
   
    _loops_file_path=loops_file_path.format(resolution=resolution)
    print(_loops_file_path)
    loops_df=pd.read_csv(_loops_file_path, sep="\t")
    loops_df["loop_number"]=np.arange(loops_df.shape[0])
    print(loops_df)
    print(loops_df.shape)

    if left_mark!=right_mark:
        print("Left and Right marks are not equal")
        _out_path=f"{out_path}/{left_mark}_{right_mark}_bidirectional"
        call(["mkdir", "-p", _out_path])
        
         # path to output file
        out_file_path_1 = f"{_out_path}/loops_left_{left_mark}_right_{right_mark}_{resolution}.bedpe"
        # path to output file
        out_file_path_2 = f"{_out_path}/loops_left_{right_mark}_right_{left_mark}_{resolution}.bedpe"        
        out_file_path_concat = f"{_out_path}/loops_{left_mark}_{right_mark}_{resolution}_bidirectional.bedpe"
        
        loops_bed_1_on_left_bed_2_on_right = intersect(loops_df, bed_1_df, bed_2_df)
        loops_bed_2_on_left_bed_1_on_right = intersect(loops_df, bed_2_df, bed_1_df)

        loops_intersect_concat = pd.concat([loops_bed_1_on_left_bed_2_on_right,loops_bed_2_on_left_bed_1_on_right])
        
        loops_bed_1_on_left_bed_2_on_right.to_csv(out_file_path_1, sep="\t", header=None, index=False)
        loops_bed_2_on_left_bed_1_on_right.to_csv(out_file_path_2, sep="\t", header=None, index=False)
        loops_intersect_concat.to_csv(out_file_path_concat, sep="\t", header=None, index=False)
    else:
        print("Left and Right marks are equal")
        _out_path=f"{out_path}/{left_mark}_{right_mark}"
        call(["mkdir", "-p", _out_path])

        # path to output file
        out_file_path = f"{_out_path}/loops_left_{left_mark}_right_{right_mark}_{resolution}.bedpe"        
        
        loops_bed_1_on_left_bed_2_on_right = intersect(loops_df, bed_1_df, bed_2_df)

        
        loops_bed_1_on_left_bed_2_on_right.to_csv(out_file_path, sep="\t", header=None, index=False)

    