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


clr = cooler.Cooler(normalized_cool_file_path)
resolution = clr.binsize
print(f"Resolution = {resolution}")
print(clr.chromnames)

print(clr.bins()[:10])
print(clr.pixels()[:10])

clr_columns = clr.bins()[:10].columns

#from https://www.biostars.org/p/9561371/
if "weight" not in clr_columns:
    weights = cooler.balance_cooler(clr)
    print(weights)
    print(weights[0].shape)
    with clr.open('r+') as f:
        f["bins"].create_dataset("weight", data=weights[0], compression="gzip", compression_opts=6)

print(clr.bins()[:10])
print(clr.pixels()[:10])
