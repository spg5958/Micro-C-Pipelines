import numpy as np
import pandas as pd
import pybedtools


loops_file_path = "../output/loops_check/loops.tsv"                   # path to loops.tsv file
bed_1_file_path = "../input/bed_files/DE_Ctrl_PRDM1.0.7rpm.bed"       # Left histone mark/TF
bed_2_file_path = "../input/bed_files/DE_Ctrl_H2AK119Ub1.all.bed"     # Right histone mark/TF
out_file_path = "../output/check.bedpe"                               # path to output file

# intersection

loops_df=pd.read_csv(loops_file_path, sep="\t")

loops_df["loop_number"]=np.arange(loops_df.shape[0])

loops_left_df=loops_df[["BIN1_CHR","BIN1_START","BIN1_END","loop_number"]]

loops_right_df=loops_df[["BIN2_CHROMOSOME","BIN2_START","BIN2_END","loop_number"]]

loops_left_bed=pybedtools.BedTool.from_dataframe(loops_left_df)

loops_right_bed=pybedtools.BedTool.from_dataframe(loops_right_df)

bed_1 = pybedtools.BedTool(bed_1_file_path)

bed_2 = pybedtools.BedTool(bed_2_file_path)

loops_left_bed_1_intersect=loops_left_bed.intersect(bed_1, u=True).to_dataframe()
loops_left_bed_1_intersect=loops_left_bed_1_intersect.rename(columns={"chrom":"BIN1_CHR", "start":"BIN1_START", "end":"BIN1_END", "name":"loop_number"})

loops_right_bed_2_intersect=loops_right_bed.intersect(bed_2, u=True).to_dataframe()
loops_right_bed_2_intersect=loops_right_bed_2_intersect.rename(columns={"chrom":"BIN2_CHROMOSOME", "start":"BIN2_START", "end":"BIN2_END", "name":"loop_number"})

loops_bed_1_on_left=loops_df.merge(loops_left_bed_1_intersect, how="inner", on=["BIN1_CHR", "BIN1_START", "BIN1_END", "loop_number"])

loops_bed_2_on_right = loops_df.merge(loops_right_bed_2_intersect, how="inner", on=["BIN2_CHROMOSOME", "BIN2_START", "BIN2_END", "loop_number"])

loops_bed_1_on_left_bed_2_on_right = loops_bed_1_on_left.merge(loops_bed_2_on_right, how="inner")

loops_bed_1_on_left_bed_2_on_right.to_csv(out_file_path, sep="\t", header=None, index=False)
