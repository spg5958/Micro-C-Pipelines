import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from coolpuppy import coolpup
import cooler
import bioframe
import multiprocessing as mp
import cooltools
from tqdm import tqdm
from coolpuppy.lib import numutils


cool_file_path="normalized.cool"            # path to normalized .cool file
bedpe_file_path="../output/check.bedpe"     # path to loops .bedpe file
flank = 25_000                              # flanking distance around loops
out_file_path="../output/check.h5"          # path to output file
out_figure_path = "../output/check.png"     # path of ouput pile-up plot
SAMPLE="DE"                                 # name of the sample


clr = cooler.Cooler(cool_file_path)

resolution = clr.binsize
print("Cooler info")
print(f"Resolution = {resolution}")
print(f"Chromonames = {clr.chromnames}")
print("Cooler bins")
print(clr.bins()[:10])

print("## get bedpe ##")
bedpe = bioframe.read_table(bedpe_file_path, schema='bedpe')
print(f"Bedpe file name = {bedpe_file_path}")
print(bedpe)

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

print("## Pileup analysis ##")
pup_bed = coolpup.pileup(clr,
                         bedpe,
                         features_format='bedpe',
                         local=False,
                         flank=flank,
                         clr_weight_name="weight",
                         expected_df=expected,
                         nproc=4
                        )

print(pup_bed)

print("## Save pileup output ##")
pup_bed.to_hdf(out_file_path, key='df', mode='w')


# PLOT

pup_bedpe = pd.read_hdf(out_file_path)
print(pup_bedpe)

insulation_strength = numutils.get_insulation_strength(pup_bedpe.loc[0, 'data'])
insulation_strength = round(insulation_strength, 2)

flank = pup_bedpe["flank"][0]
resolution = pup_bedpe["resolution"][0]

f, ax = plt.subplots()

plt.imshow(
    np.log2(pup_bedpe.loc[0, 'data']),
    vmax = 0.2,
    vmin = -0.2,
    cmap='coolwarm',
    interpolation='none')

plt.colorbar(label = 'log2 mean obs/exp')
ticks_pixels = np.linspace(0, flank*2//resolution, 5)
ticks_kbp = ((ticks_pixels-ticks_pixels[-1]/2)*resolution//1000).astype(int)
plt.xticks(ticks_pixels, ticks_kbp)
plt.yticks(ticks_pixels, ticks_kbp)
plt.xlabel('relative position, kbp')
plt.ylabel('relative position, kbp')
plt.title(SAMPLE, pad=20)
plt.text(.01, .99, str(insulation_strength), ha='left', va='top', transform=ax.transAxes, size=15)
plt.savefig(out_figure_path)