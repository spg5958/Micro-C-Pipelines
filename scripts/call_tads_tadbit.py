# TADbit normalization - https://github.com/3DGenomes/TADbit/issues/66 

import pandas as pd
from pytadbit.hic_data import HiC_data
from pytadbit.parsers.genome_parser import get_gc_content, parse_fasta
from pytadbit import Chromosome
import pytadbit
import os
print(pytadbit.__file__)
from pytadbit.hic_data import HiC_data


sample_name="DE"         # name of the sample
in_path=""               # full path to raw directory inside HiC-Pro output directory 
out_path=""              # path to output dir to tad bed files
chr_sizes_file_path=     # path to genome sizes file
chr_list=[]              # list of chromosomes
resolution_list=[]       # list of resolutions (should be present in the HiC-Pro output directory

def load_hic_data_from_hicpro(infile, resolution, chromsizes):
        print(f"Loading matrix from = {infile}")
        print(f"Loading bed file from = {chromsizes}")
        #chromosome stuff
        from collections import OrderedDict
        genome_seq = OrderedDict() #we store here chromosome bins
        #genome_seq = {} #we store here chromosome bins
        size = 0 #total size of genome in terms of bins

        sections = []
        substract = 0
        dict_sec = dict()
        #reading bed
        f = open(chromsizes, "r")
        for line in f:
                sline = line.strip("\n").split("\t")
                if sline[0] not in genome_seq:
                    substract = int(sline[3])
                binInChr = int(sline[3]) - substract + 1
                genome_seq[sline[0]] = binInChr
                #dict_sec[(sline[0], binInChr)] = int(sline[3])
                sections.append((sline[0], binInChr-1))
                size += 1
        f.close()
        dict_sec = dict([(j, i) for i, j in enumerate(sections)])
        #reading matrix
        imx = HiC_data((), size, genome_seq, dict_sec, resolution=resolution)
        f =  open(infile, "r")
        for line in f:
                sline = line.strip("\n").split("\t")
                imx[int(sline[0])-1, int(sline[1])-1] += int(sline[2])
                imx[int(sline[1])-1, int(sline[0])-1] += int(sline[2])
        f.close()
        return imx, sections, size

def get_tads(chrname, resolution, matrix_file, bed_file, sample):
    
    print(f"Chromosome = {chrname}, Resolution = {resolution}")
    
    # call TADs
    
    # load hic data into TADbit again
    hic_data, sections, size = load_hic_data_from_hicpro(matrix_file, 
                                                         resolution,
                                                         bed_file
                                                        )
    print(size)
    crm = Chromosome(chrname)
    exp_name = f"{sample}_{resolution}"
    
    crm.add_experiment(exp_name,
                       hic_data=[hic_data.get_matrix(focus=chrname)],
                       resolution=int(resolution),
                       replace = True
                      )
    
    crm.experiments[exp_name].normalize_hic()  # normalize matrix before finding TADs
    crm.find_tad([exp_name], n_cpus=4)
    return crm.experiments[exp_name].tads

def get_bed(tads_dict,_chr,resolution,chr_sizes_file_path):
    chr_sizes_df=pd.read_csv(chr_sizes_file_path, sep="\t", header=None)
    print(chr_sizes_df)
    chr_size=int(chr_sizes_df[chr_sizes_df[0]==_chr][1][0])
    print(f"chr_size = {chr_size}")
    
    bed_dict={"tad_no":[],"chr":[],"start":[],"end":[],"score":[]}
    for key,value in tads_dict.items():
        bed_dict["tad_no"].append(key)
        bed_dict["chr"].append(_chr)
        bed_dict["start"].append(value["start"]*resolution)
        bed_dict["end"].append((value["end"]+1)*resolution)
        print(value["end"])
        bed_dict["score"].append(value["score"])
    bed_df=pd.DataFrame.from_dict(bed_dict)    
    bed_df=bed_df.sort_values("tad_no")
    print(bed_df.dtypes)
   
    bed_df["start"]=bed_df["start"].astype(int)
    bed_df["end"]=bed_df["end"].astype(int)
    bed_df["end"]=bed_df["end"].clip(upper=chr_size)
    return bed_df

os.makedirs(out_path, exist_ok=True)

for CHR in chr_list:
    for RESOLUTION in resolution_list:
        os.makedirs(f"{out_path}/{CHR}/{RESOLUTION}", exist_ok=True)
        matrix_file = f"{in_path}/{RESOLUTION}/DE_merged_bio_rep_{RESOLUTION}.matrix"
        bed_file = f"{in_path}/{RESOLUTION}/DE_merged_bio_rep_{RESOLUTION}_abs.bed"
        tads_dict=get_tads(CHR, RESOLUTION, matrix_file, bed_file, f"{sample_name}")
        print(tads_dict)
        bed_df=get_bed(tads_dict, CHR, RESOLUTION, chr_sizes_file_path)
        print(bed_df)
        bed_df[["chr","start","end","score"]].to_csv(f"{out_path}/{CHR}/{RESOLUTION}/tad_domains.bed",sep="\t",header=False,index=False)
