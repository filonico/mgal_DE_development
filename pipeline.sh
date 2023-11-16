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
