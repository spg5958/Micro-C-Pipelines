import numpy as np
import pandas as pd
import pybedtools
from subprocess import call


# path to ChIP-seq .bed file
chip_seq_path = "../input/dEN_SOX17.1rpm.bed"
 
# path to chromHMM domain .bed file
chromHMM_file_path="../input/cell1_8_dense_ID.bed"
    
# chromHMM category
chromHMM_category = "Promoters"

# path to loops.tsv file
loops_file_path = "../input/loops_merged_mutli_res_Promoters_Enhancers.tsv"   

# output path
out_path = "../output/script_2"     

# path to output file
out_file_prefix = "chip-chromHMM-loop_subtypes_interesect"        

call(["mkdir", "-p", out_path])

print("chromHMM_bed_df")
chromHMM_bed_df = pd.read_csv(chromHMM_file_path, sep="\t")
print(chromHMM_bed_df)
print(chromHMM_bed_df.shape)
print(chromHMM_bed_df.columns)
print(chromHMM_bed_df["ChromHMM_Category"].unique())
chromHMM_bed_df = chromHMM_bed_df[chromHMM_bed_df["ChromHMM_Category"]==chromHMM_category]
chromHMM_bed_df = chromHMM_bed_df.iloc[:,:6]
print(chromHMM_bed_df)

print("chip_seq_bed_df")
chip_seq_bed_df = pd.read_csv(chip_seq_path, sep="\t", header=None)
print(chip_seq_bed_df)
print(chip_seq_bed_df.shape)
print(chip_seq_bed_df.columns)

print("loops_df")
print(loops_file_path)
loops_df=pd.read_csv(loops_file_path, sep="\t")
loops_df["loop_number"]=np.arange(loops_df.shape[0])
print(loops_df)
print(loops_df.shape)

# intersect ChIP-seq and chromHMM .bed file
print("bed_1_df")
chromHMM_bed = pybedtools.BedTool.from_dataframe(chromHMM_bed_df)
chip_seq_bed = pybedtools.BedTool.from_dataframe(chip_seq_bed_df)
bed_1_df = chromHMM_bed.intersect(chip_seq_bed, u=True).to_dataframe()
print(bed_1_df)

def intersect(loops_df, bed_1_df, bed_2_df):
    
    # intersection
    loops_left_df=loops_df[["chrom1","start1","end1","loop_number"]]

    loops_right_df=loops_df[["chrom2","start2","end2","loop_number"]]

    loops_left_bed=pybedtools.BedTool.from_dataframe(loops_left_df)

    loops_right_bed=pybedtools.BedTool.from_dataframe(loops_right_df)

    # https://github.com/daler/pybedtools/issues/238
    bed_1 = pybedtools.BedTool.from_dataframe(bed_1_df)
    bed_2 = pybedtools.BedTool.from_dataframe(bed_2_df)

    loops_left_bed_1_intersect=loops_left_bed.intersect(bed_1, u=True).to_dataframe()
    loops_left_bed_1_intersect=loops_left_bed_1_intersect.rename(columns={"chrom":"chrom1", "start":"start1", "end":"end1", "name":"loop_number"})

    loops_right_bed_2_intersect=loops_right_bed.intersect(bed_2, u=True).to_dataframe()
    loops_right_bed_2_intersect=loops_right_bed_2_intersect.rename(columns={"chrom":"chrom2", "start":"start2", "end":"end2", "name":"loop_number"})

    loops_bed_1_on_left=loops_df.merge(loops_left_bed_1_intersect, how="inner", on=["chrom1", "start1", "end1", "loop_number"])

    loops_bed_2_on_right = loops_df.merge(loops_right_bed_2_intersect, how="inner", on=["chrom2", "start2", "end2", "loop_number"])

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
                
    loops_bed_1_on_left.to_csv(f"{out_path}/{out_file_prefix}_left.bedpe", sep="\t", header=None, index=False)
    loops_bed_2_on_right.to_csv(f"{out_path}/{out_file_prefix}_right.bedpe", sep="\t", header=None, index=False)
    loops_bed_1_on_left_bed_2_on_right.to_csv(f"{out_path}/{out_file_prefix}_left_right.bedpe", sep="\t", header=None, index=False)


intersect(loops_df, bed_1_df, bed_1_df)
