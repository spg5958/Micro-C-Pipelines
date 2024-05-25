# [DOCUMENTATION & SOME SCRIPTS ARE UNDER DEVELOPMENT]

# Pipeline 1

## Step 1 - Process Individual Replicates
- Arrange the FASTQ files according to the following directory structure. HiC-Pro considers all readsr within one input folder as one sample
```
|--PATH_TO_FASTQ_FILES
   |--sample1
      |--file1_R1.fastq.gz
      |--file1_R2.fastq.gz
      |--file2_R1.fastq.gz
      |--file2_R2.fastq.gz
      |-- ...
   |--sample2
     |-- file1_R1.fastq.gz
     |-- file1_R2.fastq.gz
```
- Generate annotation files:
  - chromosomes sizes file: Two-column tab-separated text file containing chromosome names and sizes. 
e.g.: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes 
   - The bowtie2 indexes: See the bowtie2 manual page for details about how to create such indexes. 
   - (HiC only) A BED file of the restriction fragments after digestion: This file depends on both the restriction enzyme and the reference genome. See the HiC-Pro annotation for details about how to generate this file.
 
- Setup the configuration file: <br>
  Copy and edit the 'config-hicpro.txt' file in your local directory. Modify the options in the configuration file as required. Below are some important options to consider; for more detailed explanations of these options, refer to the HiC-Pro manual. 
  - PAIR1_EXT = Keyword for first mate detection. Default: _R1 
  - PAIR2_EXT = Keyword for second mate detection. Default: _R2 
  - REFERENCE_GENOME = Reference genome prefix used for genome indexes. Default: hg19 
  - GENOME_SIZE = Full path to chromosome size file. 
  - MIN_CIS_DIST4 = [Important option for Micro-C data] Filter short range contact below the specified distance. Example: 1000 
  - BIN_SIZE = Resolution of contact maps to generate (space separated). Default: 20000 40000 150000 500000 1000000 
  - (Hi-C only) GENOME_FRAGMENT = Full path to BED file with restriction fragments. 
  - (Hi-C only) LIGATION_SITE = Ligation site sequence(s) used for reads trimming.      - Depends on the fill in strategy. Note that multiple ligation sites can be specified (comma separated). Example: AAGCTAGCTT 
  - (Hi-C only) MIN_FRAG_SIZE = Maximum size of restriction fragments to consider for the Hi-C processing. Example: 100 
  - (Hi-C only) MAX_FRAG_SIZE = Maximum size of restriction fragments to consider for the Hi-C processing. Example: 100000 
  - (Hi-C only) MIN_INSERT_SIZE = Minimum sequenced insert size. Shorter 3C products are discarded. Example: 100 
  - (Hi-C only) MAX_INSERT_SIZE = Maximum sequenced insert size. Larger 3C products are discarded. Example: 600

- Run HiC-pro.
```
Run:
MY_INSTALL_PATH/bin/HiC-Pro -i FULL_PATH_TO_DATA_FOLDER -o FULL_PATH_TO_OUTPUTS -c MY_LOCAL_CONFIG_FILE
```
## Step 2 - Merge Replicates
