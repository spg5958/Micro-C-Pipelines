import pandas as pd
from subprocess import call
import os
import subprocess

# resoution and corresponding .cool files location
res_dict = {
        "500000":"../output/DE_merged_bio_rep_500000_normalized.cool",
        "600000":"../output/DE_merged_bio_rep_600000_normalized.cool",
        "700000":"../output/DE_merged_bio_rep_700000_normalized.cool",
        "800000":"../output/DE_merged_bio_rep_800000_normalized.cool"
        }

# chromosome list
chr_list = "chr1"

# path to output directory
out_path="/storage/home/spg5958/group/lab/siddharth/projects/micro_c_analysis_Iwafuchi_lab/MCAI64-R1-call_loops_at_multiple_res_and_combine/output/"

# path to mustache.py
mustache_path = "/storage/home/spg5958/group/lab/siddharth/software/mustache/mustache/mustache.py"

if chr_list=="all":
    for res, cool_file_path in res_dict.items():
        call([mustache_path, "-f", cool_file_path, "-r", res, "-st", "0.7", "-pt", "0.5", "-p", "4", "-o", f"{out_path}/loops_chr_all_{res}.tsv"])
else:
    for res, cool_file_path in res_dict.items():
        output=subprocess.run([f"{mustache_path} -f {cool_file_path} -ch {chr_list} -r {res} -st 0.7 -pt 0.5 -p 4 -o {out_path}/loops_{chr_list.replace(' ','_')}_{res}.tsv"], shell=True, capture_output=True, text=True)
        print(output.stdout)
        print(output.stderr)

for file_name in os.listdir(out_path):
    if file_name.endswith(".tsv"):
        print(file_name)
        df=pd.read_csv(f"{out_path}/{file_name}",sep="\t")
        out_df=df.iloc[:,:6]
        print(out_df)
        out_df.to_csv(f"{out_path}/{file_name.split('.')[0]}.bedpe", sep="\t", header=None, index=False)

