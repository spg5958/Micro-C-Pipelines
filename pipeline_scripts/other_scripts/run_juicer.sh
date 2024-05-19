#!/bin/bash
#SBATCH --account=sam77_h
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=160GB
#SBATCH --partition=sla-prio
#SBATCH --job-name=juicer
#SBATCH --output=slurm-%x-%j.out
umask 007

#expects $data

cd /storage/group/sam77/default/seqalign/reads/Iwafuchi/20221209_MicroC/Seq018_MicroC_FASTQ/juicer_analysis

./scripts/juicer.sh -g hg38 -d "/storage/group/sam77/default/seqalign/reads/Iwafuchi/20221209_MicroC/Seq018_MicroC_FASTQ/juicer_analysis/"$data -D /storage/group/sam77/default/seqalign/reads/Iwafuchi/20221209_MicroC/Seq018_MicroC_FASTQ/juicer_analysis/ -t 10  -z references/hg38/hg38.fa -p references/hg38/chrom.sizes
