#!/bin/bash

##############################
#     Download raw reads     #
##############################

mkdir -p 01_raw_reads/01_fastqc

# download reads and perform quality check
# REQUIRES: conda_env/ncbi_env.yml
python3 scripts/01_download_reads.py -i 00_input/SRA_toDownload.ls -o 01_raw_reads

# aggregate fastqc report into a single html file
# REQUIRES: conda_envs/multiqc_env.yml
multiqc -o 01_raw_reads/01_fastqc/ 01_raw_reads/01_fastqc/


######################
#     TRIM READS     #
######################

# create a directory to store trimmed reads and quality control results
mkdir -p 02_trimmed_reads/01_fastqc

# trim reads
# REQUIRES: conda_envs/readtrimming_env.yml
for i in 01_raw_reads/SRR*; do python3 scripts/02_trim_reads.py -d $i -adapt 00_input_files/contaminants2trimm.fa; done

# aggregate fastqc report into a single html file
multiqc -o 02_trimmed_reads/01_fastqc 02_trimmed_reads/01_fastqc


###############################
#     MAP READS WITH STAR     #
###############################

# index the reference genome and map reads using STAR
# REQUIRES: conda_envs/mapreads_star_env.yml
bash scripts/03a_map_reads_star.sh

# get read mapping raw counts
# REQUIRES: conda_envs/countreads_stringtie_env.yml
bash scripts/05a_getcounts_stringtie.sh


#################################
#     MAP READS WITH BOWTIE     #
#################################

# index transcriptome, map reads, convert sam to bam and get raw counts statistics
# REQUIRES: conda_envs/mapreads_bowtie_env.yml
bash scripts/03b_map_reads_bowtie.sh

# merge raw count statistics in one file
scripts/05b_merge_rawmappings_bowtie.sh

# perform PCA analysis on read counts
# REQUIRES: conda_envs/R_env.yml
mkdir 04_PCA_readcounts
Rscript scripts/06_plotPCA_readcounts.R


###########################
#     RUN DE ANALYSIS     #
###########################

mkdir 05_masigpro_analysis

# normalize read counts and run a PCA analysis
Rscript scripts/07_DE_timeseries.R
