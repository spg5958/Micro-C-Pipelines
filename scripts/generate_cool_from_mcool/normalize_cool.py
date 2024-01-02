import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cooler
import multiprocessing as mp
import cooltools
from tqdm import tqdm
from subprocess import call
import sys

normalized_cool_file_path=sys.argv[1]
print(normalized_cool_file_path)

#"../output/cooler_balance_python/DE_merged_bio_rep_500000_normalized.cool"

clr = cooler.Cooler(normalized_cool_file_path)
resolution = clr.binsize
print(f"Resolution = {resolution}")
print(clr.chromnames)

print(clr.bins()[:10])
print(clr.pixels()[:10])

clr_columns = clr.bins()[:10].columns

if "weight" not in clr_columns:
    pool = mp.Pool(4)
    weights = cooler.balance_cooler(clr,map=pool.map)
    print(weights)
    print(weights[0].shape)
    with clr.open('r+') as f:
        f["bins"].create_dataset("weight", data=weights[0], compression="gzip", compression_opts=6)

print(clr.bins()[:10])
print(clr.pixels()[:10])
