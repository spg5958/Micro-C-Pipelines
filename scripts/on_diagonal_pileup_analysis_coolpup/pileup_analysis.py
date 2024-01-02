import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from coolpuppy import coolpup
import cooler
import bioframe
import multiprocessing as mp
import cooltools
from tqdm import tqdm


cool_file_path="../../input/DE_merged_bio_rep_1000.cool"
bed_file_name="PRDM1_vs_H2AK119Ub1-H2AK119Ub1.no_PRDM1.bed"
flank = 25_000
out_file_path="../../output/PRDM1_vs_H2AK119Ub1-H2AK119Ub1.no_PRDM1/pileup_analysis_out_PRDM1_vs_H2AK119Ub1-H2AK119Ub1.no_PRDM1.h5"

clr = cooler.Cooler(cool_file_path)

resolution = clr.binsize
print("Cooler info")
print(f"Resolution = {resolution}")
print(f"Chromonames = {clr.chromnames}")
print("Cooler bins")
print(clr.bins()[:10])

clr_columns = clr.bins()[:10].columns

if "weight" not in clr_columns:
    pool = mp.Pool(4)
    weights = cooler.balance_cooler(clr,map=pool.map)
    print(weights)
    print(weights[0].shape)
    with clr.open('r+') as f:
        f["bins"].create_dataset("weight", data=weights[0], compression="gzip", compression_opts=6)
        
print(clr.bins()[:10])

print("## Calcualtion expected interactions ##")
expected = cooltools.expected_cis(clr, nproc=4, chunksize=1_000_000)
print(expected)

print("## get bed ##")
bed = bioframe.read_table(f"../../input/bed_files/{bed_file_name}", schema='bed')
print(f"Bed file name = {bed_file_name}")
print(bed)

print("## Pileup analysis ##")
pup_bed = coolpup.pileup(clr,
                         bed,
                         features_format='bed',
                         local=True,
                         flank=flank,
                         clr_weight_name="weight",
                         expected_df=expected,
                         nproc=4
                        )

print(pup_bed)

print("## Save pileup output ##")
pup_bed.to_hdf(out_file_path, key='df', mode='w')
