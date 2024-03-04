import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from coolpuppy.lib import numutils
import os


h5_files_dir_path="../output"
out_figure_path = "../output"

file_list=[]
for file in os.listdir(h5_files_dir_path):
    if file.endswith('.h5'):
        print(file)
        file_list.append(file)

for file in file_list:
    _file_name = file.split(".")[0]
    pup_bedpe = pd.read_hdf(f"{h5_files_dir_path}/{file}")
    print(pup_bedpe)

    insulation_strength = numutils.get_insulation_strength(pup_bedpe.loc[0, 'data'])
    insulation_strength = round(insulation_strength, 2)

    flank = pup_bedpe["flank"][0]
    resolution = pup_bedpe["resolution"][0]

    f, ax = plt.subplots()

    plt.imshow(
        np.log2(pup_bedpe.loc[0, 'data']),
    #     vmax = 0.2,
    #     vmin = -0.2,
        cmap='coolwarm',
        interpolation='none')

    plt.colorbar(label = 'log2 mean obs/exp')
    ticks_pixels = np.linspace(0, flank*2//resolution, 5)
    ticks_kbp = ((ticks_pixels-ticks_pixels[-1]/2)*resolution//1000).astype(int)
    plt.xticks(ticks_pixels, ticks_kbp)
    plt.yticks(ticks_pixels, ticks_kbp)
    plt.xlabel('relative position, kbp')
    plt.ylabel('relative position, kbp')
    plt.title(_file_name, pad=20)
    plt.text(.01, .99, str(insulation_strength), ha='left', va='top', transform=ax.transAxes, size=15)
    plt.savefig(f"{out_figure_path}/{_file_name}.png")
