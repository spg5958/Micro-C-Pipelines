import pandas as pd
import numpy as np
import bioframe
import scipy

pca_bedgraph_file_path="../input/MCAI34-R2-output/pca1_hFG-lib548+lib549_merged_chr1_1000000.bedgraph"
hg38_fasta_file_path="/storage/home/spg5958/group/genomes/hg38/hg38.analysisSet.fa"
out_path = "../output"

df = pd.read_table(pca_bedgraph_file_path, header=None)
fasta_records=bioframe.load_fasta(hg38_fasta_file_path)

print(df)

# for one chromosome

compartments_bed_dict={}
compartment_idx = 0
prev_val = None
for i,curr_row in df.iterrows():
    curr_val = curr_row[3]
    if len(compartments_bed_dict) == 0:
        compartments_bed_dict[compartment_idx] = {"chr":curr_row[0] , "start":curr_row[1], "end":curr_row[2], "pca_sign":np.sign(curr_val)}
    elif prev_val*curr_val >= 0:
        compartments_bed_dict[compartment_idx]["end"] = curr_row[2]
    elif prev_val*curr_val < 0:
        compartment_idx += 1
        compartments_bed_dict[compartment_idx] = {"chr":curr_row[0] , "start":curr_row[1], "end":curr_row[2], "pca_sign":np.sign(curr_val)}
    prev_val = curr_val
    
compartments_bed_df = pd.DataFrame()
for i in range(len(compartments_bed_dict)):
    compartments_bed_df = compartments_bed_df.append(compartments_bed_dict[i], ignore_index=True)
print(compartments_bed_df)   


df1=df.rename(columns={0:"chrom", 1:"start", 2:"end", 3:"pca"})
df1=bioframe.frac_gc(df1,fasta_records)
print(df1)
print(df)

print(df1["GC"].isna().sum())
df1=df1.dropna(subset=["GC"])
r, p = scipy.stats.spearmanr(df1["pca"],df1["GC"])
print(r)

compartment_list=[]
if r>0:
    for i in compartments_bed_df["pca_sign"]:
        if i>=0:
            compartment_list.append("A")
        else:
            compartment_list.append("B")       
else:
    for i in compartments_bed_df["pca_sign"]:
        if i>=0:
            compartment_list.append("B")
        else:
            compartment_list.append("A")     
compartments_bed_df["compartment"]=compartment_list

compartments_bed_df.to_csv(f"{out_path}/compartments_bed_df.csv",index=False)
print(compartments_bed_df)

