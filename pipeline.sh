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


#####################
#     MAP READS     #
#####################

# Genome assembly and annotation were downloaded from https://www.ncbi.nlm.nih.gov/assembly/GCA_900618805.1
# The GFF genome annotation file was converted into a GTF file by running the following command:
# agat_convert_sp_gff2gtf.pl --gff GCA.900618805.1_mgal_genomic.gff.gz -o GCA.900618805.1_mgal_genomic.gtf

# index the reference genome and map reads using STAR
# REQUIRES: conda_envs/mapreads_star_env.yml
bash scripts/03_map_reads_star.sh

# get read mapping raw counts
# REQUIRES: conda_envs/countgreads_stringtie_env.yml
bash scripts/05_getcounts_stringtie.sh
