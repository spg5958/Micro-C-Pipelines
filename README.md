# [DOCUMENTATION & SOME SCRIPTS ARE UNDER DEVELOPMENT]

<br>
<br>

# Pipeline 1 - Pipeline for generating contact matrix from FASTQ files using HiC-Pro

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
To merge biological replicates after merging technical replicates, we can rerun HiC-Pro starting from the validPairs files generated during the previous analysis. There's no need to start from the FASTQ files again.

- Arrange validPairs files according to the following directory structure.
  ```
  |--PATH_TO_validPairs_FILES
     |--validPairs
        |--validPairs tech_sample1
        |--validPairs tech_sample2
  ```
- Rerun Hic-Pro in stepwise mode.
  ```
  MY_INSTALL_PATH/bin/HiC-Pro -i FULL_PATH_TO_VALIDPAIRS_DATA -o FULL_PATH_TO_OUTPUTS -c MY_LOCAL_CONFIG_FILE -s merge_persample -s build_contact_maps -s ice_norm
  ```
  
## Step 3 - Generate .cool Files
   HiC-Pro outputs contact matrices in .matrix format. However, .cool files are generally accepted by most HiC analysis software and are also easier to handle than other formats.
   Therefore, I convert .matrix files to .cool files using the following 'hicpro2higlass.sh' script that comes with HiC-Pro. 
   ```
   Run:
      HICPRO_PATH/bin/utils/hicpro2higlass.sh -i MATRIX_FILE -r RESOLUTION -c CHROMSIZES_FILE -o OUTPUT_PATH -p NUM_PROC
   ```

## Step 4 - Normalize .cool Files
   Although HiC-Pro can generate normalized .cool files, I prefer to normalize the raw .cool files separately. I use `pipeline_scripts/4_normalize_cool/normalize_cool.py` Python script for normalizing the .cool files. This script uses Iterative correction matrix balancing method from cooler package. The full path to the .cool file should be provided at line number before running.

<br>
<br>

# Pipeline 2 - Pipeline for Identification of A/B Compartments 
This pipeline is based on the protocol outlined by Miura et al. in their book1. It utilizes HiCExplorer commands as detailed in the documentation2 and generates A/B compartments starting from a raw .cool matrix. While input for pipeline is a raw .cool matrix, but you can start from a normalized .cool matrix and bypass the step 1 below.

1) Normalize .cool file:
   - (optional) Normalize to equal level of read coverage or value ranges: If you want to compare different Hi-C interaction matrices (samples/replicates), the matrices need to be normalized to equal level of read coverage or value ranges. This can be achieved with ‘hicNormalize’ command as follows. Please refer to hicNormalize documentation for further details.
     ```
     Run:
       hicNormalize -m matrix.cool --normalize smallest -o normalized_matrix.cool
     ```
   - Correct HiC matrix using Iterative Correction (ICE): The Hi-C matrix has to be corrected to remove GC, open chromatin biases and, most importantly, to normalize the number of restriction sites per bin4. However, for this method to work correctly, bins with zero reads assigned to them should be removed as they cannot be corrected. Also, bins with the low number of reads should be removed. Bins with an extremely high number of reads can also be removed from the correction as they may represent copy number variations. To aid in the identification of bins (threshold values) with low and high read coverage, the histogram (diagnostic plot) of the number of reads can be plotted together with the Median Absolute Deviation (MAD) using the ‘hicCorrectMatrix’5 command as shown below.
     ```
     Run:
       hicCorrectMatrix diagnostic_plot –m matrix.cool -o plot_file.png --chromosomes chr
     ```
     Once the thresholds have been decided, the matrix can be corrected as follows.
     ```
     Run:
       hicCorrectMatrix correct -m matrix.cool --correctionMethod ICE –chromosomes chr_list --iterNum 5000 -o corrected.cool --filterThreshold -1.5 5
     ```
  2) Generate observed-over-expected matrix:
     To compute the A/B compartments the matrix needs to be transformed to an observed-over-expected matrix using. ‘hicTransform’7 command.
     ```
     Run:
       hicTransform -m corrected.cool --method obs_exp -o out.cool --chromosomes chr_list
     ```
  3) Calculate Pearson correlation matrix:
     Next, we transform the observed-over-expected matrix into a Pearson correlation matrix by calculating the Pearson correlation values between all pairs of rows and columns as follows. 
     ```
     Run:
       hicTransform -m obs_exp.cool --method pearson -o pearson.cool --chromosomes chr_list 
     ```
  4) Compute first principal component of the Pearson correlation matrix:
     The first principal component (PC) of the Pearson correlation matrix has been shown to correlate with A/B compartments. The first PC track (in bedgraph format) is computed using ‘hicPCA’ command from HiCExplorer as shown below.
     ```
     Run:
       hicPCA -m pearson.cool -o pca1.bedgraph --format bedgraph -we 1 --chromosomes chr_list 
     ```
  5) Calculate gene density:
     Unfortunately, PCA does not provide information about the relationship between the sign of the first principal component (PC) and the A/B compartments. Therefore, it is recommended to use external information, such as histone marks, to identify A/B compartments. One of the straightforward methods to determine A/B compartments is by calculating gene density. The A compartment (active) tends to have a high gene density. It is estimated that regions with high GC content have a higher relative gene density compared to regions with lower GC content10,11. The following Python script calculates the GC content (i.e., gene density) track for given chromosomes. This script requires the TADbit package. Please provide appropriate input arguments at the beginning of the script (lines 9-11).
     ```
     # conda env = tadbit 

    import numpy as np 
    import pandas as pd 
    from pytadbit.parsers.genome_parser import get_gc_content, parse_fasta 
    import os.path 
    from subprocess import call 

    hic_pro_abs_bed_file=""   # full path to abs.bed file generated by HiC-Pro 
    genome_fasta_file=""      # full path to genome fasta file 
    out_bedgraph_file=""      # full path to output bedgraph file  
        
    df = pd.read_csv(hic_pro_abs_bed_file,sep="\t",header=None) 
    df_chr=df[df[0]==CHR] 
    print(df_chr) 

    rich_in_A = get_gc_content(parse_fasta(genome_fasta_file),chromosomes=[CHR],resolution=RESOLUTION, n_cpus=4) 
    print(rich_in_A) 
    print(len(rich_in_A)) 
    df_chr[4]=rich_in_A 
    df_chr=df_chr.dropna(subset=4) 
    df_out=df_chr[[0,1,2,4]] 
    print(df_out) 
    
    df_out.to_csv(out_bedgraph_file, header=False,index=False,sep="\t") 
    ```
  6) Plot A/B compartments:
     A/B compartments can be visualized using ‘pyGenomeTracks’13 package. To begin, you need to create a track file that includes the first principal component (PC), gene density, and the contact matrix. Below is an example of how to structure your track file, along with the command to plot the A/B compartments using this file.
     Example track.ini file: 
  ```
  [x-axis] 
  where = top 

  [spacer] 
  
  [pca1 bedgraph] 
  file = pca1_Seq021_MicroC_merged_rep_chr1_500000.bedgraph 
  color = blue 
  height = 5 
  transform = no 
  negative_color = red 
  height = 2 
  title = PCA 1 

  [rich_in_A.bedgraph] 
  file = rich_in_A_Seq021_MicroC_merged_rep_chr1_500000.bedgraph 
  color = orange 
  type = line:2 
  height = 2 
  title = GC content 

  [spacer] 

  [hic matrix] 
  file = obs_exp_Seq021_MicroC_merged_rep_chr1_500000.cool 
  depth =  250_000_000 
  file_type = hic_matrix_square 
  colormap = coolwarm 
  transform = log 
  min_value = -2 
  max_value = 2 
  show_masked_bins = true 
  title = log(O/E)
  ``` 
